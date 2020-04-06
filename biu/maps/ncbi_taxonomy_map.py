from ..structures import Dataset2
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule('pandas')

class NCBITaxonomyMap(Dataset2):
    """
    A mapping object for NCBI Taxonomy items.
    It is quite large so it might take some time to load...
    
    Example usage:
    --------------
    tm = biu.maps.TaxMap(grch37=True)
    
    mapping = tm.tax_id[9606]
    mapping.scientific_name # Homo sapiens
    
    tm.scientific_name['Homo sapiens'].tax_id # -> 9606
    """
    def __init__(self, *pargs, **kwargs):
        super(NCBITaxonomyMap, self).__init__("NCBITaxonomyMap", *pargs, **kwargs)

        interesting_fields    = ['tax_id', 'scientific name', 'common name', 'synonym', 'genbank common name' ]
        interesting_fields_ul = ['tax_id', 'scientific_name', 'common_name', 'synonym', 'genbank_common_name' ]
        index_fields          = ['tax_id', 'scientific_name', 'common_name', 'synonym' ]
        
        def pivot(infile, outfile):
            D = pd.read_csv(infile[0], names=['tax_id', 'name_txt', 'unique_name', 'name_class'], sep='\t')
            D = D[tax.names.name_class.isin(interesting_fields)]\
                    .drop_duplicates(['tax_id', 'name_class'])\
                    .pivot(index='tax_id',
                           columns='name_class',
                           values='name_txt')\
                    .reset_index()\
                    .rename(columns={c: c.replace(' ', '_') for c in interesting_fields})[index_fields]
            D.to_pickle(outfile)
            
            return biu.utils.Acquire2.STATUS_SUCCESS
        #edef
        
        F = utils.Acquire2()\
            .curl('ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz')\
            .gunzip()\
            .untar()\
            .select('names.dmp')\
            .cmd("sed {} -e 's/\t|\t/\t/g' -e 's/\t|$//'", '{}')\
            .func(pivot)
        

        self._obj.add_file('table.pkl', F)
        self._obj.register('_table', ['table.pkl'], lambda f: pd.read_pickle(f['table.pkl']))

        for idx, col in enumerate(interesting_fields_ul):
            if col in ['tax_id', 'scientific_name', 'synonym', 'common_name' ]:
                self._obj.register(col, [], lambda x, idx=idx: formats.MappingIndex(self._table, idx))
            #fi
        #efor
    #edef
#eclass