from ..structures import Dataset2
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule('pandas')

###############################################################################

#https://www.ensembl.org/info/data/biomart/biomart_restful.html#wget

class BioMart(Dataset2):
    """
    An interface to the BioMart identifier mapping database.
    
    Two versions are defined:
        * hsapiens_gene_trans_prot_geneid_hgnc
            Provides HUMAN ensembl geneid, transcript and protein ids, and the gene mappings to entrez, ucsc and HGNC.
        * mmusculus_gene_trans_prot_hgnc
            Provides MOUSE ensembl geneid, transcript and protein ids, and the gene mappings to HGNC.

        
    Example usage:
    --------------
    bm = biu.maps.BioMart(grch37=True)
    bm.hgnc_symbol['MTOR']
    
    bm_other = biu.maps.BioMart(database="mmusculus_gene_ensembl", attributes=["ensembl_gene_id", "hgnc_symbol"])
    bm_other = bm_other.hgnc_symbol["MTOR"]
    
    """

    versions = {
        "hsapiens_gene_trans_prot_geneid_hgnc" : {
            "database" : "hsapiens_gene_ensembl",
            "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'entrezgene', 'uniprotswissprot', 'hgnc_symbol' ]
        },
        "mmusculus_gene_trans_prot_hgnc" : {
            "database" : "mmusculus_gene_ensembl",
            "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_symbol' ]
        }
    }

    def __init__(self, version=None, database=None, attributes=None, grch37=False, *pargs, **kwargs):
        """
        Generate a BioMart Map.
        Inputs:
        -------
          version: String. Predefined dataset to use. See biu.maps.BioMart.versions()
          database: String. If no version is defined, use this database. See the XML structure of the biomart query to see the database
          attributes: List of Strings. Attributes to retrieve. Note that you can only request up to 3 external attributes
          grch37: Boolean. Use the GRCH37 build (default False)
          pargs, kwargs: See options for biu.structures.Dataset2
        Outputs:
          A Biomart Map
        """
        
        if (version is None) and (database is None) and (attributes is None):
            version = list(self.versions.keys())[0]
            database = self.versions[version]["database"]
            attributes = self.versions[version]["attributes"]
        elif (version is None) and (database is not None) and (attributes is not None):
            pass
            #version = "%s.%s%s" % (database, 'grch37.' if self.__grch37 else '', '_'.join(attributes))
        elif version is not None:
            database = self.versions[version]["database"]
            attributes = self.versions[version]["attributes"]
        else:
            raise ValueError('You must either specify a version, or a dataset and attributes.')
            return None
        #fi
        
        version = "%s.%s%s" % (database, 'grch37.' if grch37 else '', ','.join(attributes))
        
        super(BioMart, self).__init__("BioMart/%s" % version, *pargs, **kwargs)
    
        url = "http://grch37.ensembl.org/biomart/martservice?query=" if grch37 else 'http://www.ensembl.org/biomart/martservice?query='
    
        def genQuery(url, database, attributes):
            query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
            query += '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "0" datasetConfigVersion = "0.6" >'
            query += '<Dataset name = "%s" interface = "default" >' % database
            query += ''.join(['<Attribute name = "%s" />' % attr for attr in attributes ])
            query += '</Dataset></Query>'
      
            return url + query
        #edef
        
        db = utils.Acquire2().wget(genQuery(url, database, attributes))
        self._obj.add_file("biomart.tsv", db)
        
        self._obj.register('_table', ["biomart.tsv"], lambda f: pd.read_csv(f["biomart.tsv"], sep='\t',
                                                                    index_col=False, dtype=object,
                                                                    names=attributes))
        
        ##################
        # Here there is a little trick to make sure that items like this are properly registered:
        # lambda f, const1=const1, const2=const2: func(f, const1, const2)
        # Internally in the Dataset, it is called as func(f), and the consts are constant.
        ##################
        
        for idx, attr in enumerate(attributes):
            self._obj.register(attr, ["biomart.tsv"], lambda f, idx=idx: formats.MappingIndex(self._table, idx))
        #efor

    #edef
#eclass

