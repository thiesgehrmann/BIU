from ..structures import Dataset
from .. import formats
from .. import utils
from ..config import settings

###############################################################################

#https://www.ensembl.org/info/data/biomart/biomart_restful.html#wget

class BioMart(Dataset):

  versions = {
    "hsapiens_gene_trans_prot_geneid_hgnc" : {
      "database" : "hsapiens_gene_ensembl",
      "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'entrezgene', 'ucsc', 'hgnc_symbol' ]
    },
    "mmusculus_gene_trans_prot_hgnc" : {
      "database" : "hsapiens_gene_ensembl",
      "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_symbol' ]
    }

    
  }

  def __init__(self, version=None, database=None, attributes=None, where=None, grch37=False, **kwargs):
    self.__grch37 = grch37
    url = "http://grch37.ensembl.org/biomart/martservice?query=" if grch37 else 'http://www.ensembl.org/biomart/martservice?query='

    if (version is None) and (database is not None) and (attributes is not None):
      version = list(self.versions.keys())[0]
    if (version is None) and (database is not None) and (attributes is not None):
      url = url if (url is not None) else 'http://www.ensembl.org/biomart/martservice?query='
      version = "%s.%s.%s" % ([ p for p in url.split('/') if '.org' in p ][0], database, '_'.join(attributes))
    else:
      url = self.versions[version]["url"]
      database = self.versions[version]["database"]
      attributes = self.versions[version]["attributes"]
    #fi

    fileIndex = self.__genFileIndex(version, url, database, attributes, where)
    Dataset.__init__(self, fileIndex, **kwargs)

    self.__url = url
    self.__database = database
    self.__attributes = attributes
    self.__version = version

    for idx, attr in enumerate(self.__attributes):
      self._registerObject(attr, formats.TSVIndex, [ 'db' ], fileIndex['db'].path, idx, names=self.__attributes, delimiter='\t')
    #efor
  #edef

  def __genFileIndex(self, version, url, database, attributes, where):

    def genQuery():
      query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
      query += '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "0" datasetConfigVersion = "0.6" >'
      query += '<Dataset name = "%s" interface = "default" >' % database
      query += ''.join(['<Attribute name = "%s" />' % attr for attr in attributes ])
      query += '</Dataset></Query>'
  
      return url + query
    #edef

    finalPath = '%s/bioMart_%s' % ( (settings.getWhere() if where is None else where), version)
    files = {}
    files['db'] = utils.Acquire(where=where).wget(genQuery()).finalize('%s/data.tsv' % version)
    return files
  #edef

  def __genQuery(self, url, database, attributes):
    url = self.versions[version]['url']
    database = self.versions[version]['database']
    attributes = self.versions[version]['attributes']

    query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
    query += '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "0" datasetConfigVersion = "0.6" >'
    query += '<Dataset name = "%s" interface = "default" >' % database
    query += ''.join(['<Attribute name = "%s" />' % attr for attr in attributes ])
    query += '</Dataset></Query>'

    return url + query
  #edef

