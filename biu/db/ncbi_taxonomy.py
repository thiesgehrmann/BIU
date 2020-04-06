from ..structures import Dataset2
from .. import utils

pd = utils.py.loadExternalModule('pandas')


class NCBITaxonomy(Dataset2):
    """
    The NCBI Taxonomy dataset.
    Available objects:
     * nodes
     * names
     * gencode
     
    See documentation about each table at
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
    
    """
    def __init__(self, *pargs, **kwargs):
        super(NCBITaxonomy, self).__init__("NCBITaxonomy", *pargs, **kwargs)
        F = utils.Acquire2()\
                .curl('ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz')\
                .gunzip()\
                .untar()
        
        file_columns = {
            'nodes.dmp' : ['tax_id', 'parent_tax_id', 'rank', 'embl_code', 'division_id', 'inherited_div_flag', 'genetic_code_id', 'inherited_GC_flag', 'mitochondrial_genetic_code', 'inherited_MGC_flag', 'GenBank_hidden_flag', 'hidden_subtree_root_flag', 'comments' ],
            'names.dmp' : [ 'tax_id', 'name_txt', 'unique_name', 'name_class' ],
            'gencode.dmp' : [ 'genetic_code_id', 'abbreviation', 'name', 'cde', 'starts' ]
        }
        
        for file in file_columns:
            self._obj.add_file(file, F.select(file).cmd("sed {} -e 's/\t|\t/\t/g' -e 's/\t|$//'", '{}'))
            
            self._obj.register(file.split('.')[0], [file],
                               lambda x, file=file: pd.read_csv(x[file], names=file_columns[file], sep='\t'))
        #efor
    #edef
#eclass