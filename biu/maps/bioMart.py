from ..structures import Dataset
from .. import formats
from .. import utils
from ..config import settings

###############################################################################

#https://www.ensembl.org/info/data/biomart/biomart_restful.html#wget

class BioMart(Dataset):

  __versions = {
    "hsapiens_gene_trans_prot_geneid_hgnc" : {
      "database" : "hsapiens_gene_ensembl",
      "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'entrezgene', 'ucsc', 'hgnc_symbol' ]
    },
    "mmusculus_gene_trans_prot_hgnc" : {
      "database" : "mmusculus_gene_ensembl",
      "attributes" : [ 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_symbol' ]
    }

    
  }

  @classmethod
  def versions(cls):
    print("Pre-defined BioMart queries:")
    for version in cls.__versions:
      print(" * %s" % version)
      print("  - database: %s" % cls.__versions[version]['database'])
      print("  - attributes: %s" % ','.join(cls.__versions[version]['attributes']))
    #efor
  #edef

  def __init__(self, version=None, database=None, attributes=None, where=None, grch37=False, **kwargs):
    self.__grch37 = grch37
    url = "http://grch37.ensembl.org/biomart/martservice?query=" if grch37 else 'http://www.ensembl.org/biomart/martservice?query='

    if (version is None) and (database is None) and (attributes is None):
      version = list(self.__versions.keys())[0]
      database = self.__versions[version]["database"]
      attributes = self.__versions[version]["attributes"]
    elif (version is None) and (database is not None) and (attributes is not None):
      version = "%s.%s%s" % (database, 'grch37.' if self.__grch37 else '', '_'.join(attributes))
    elif version is not None:
      database = self.__versions[version]["database"]
      attributes = self.__versions[version]["attributes"]
    else:
      utils.msg.error('You must either specify a version, or a dataset and attributes.')
      return None
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

    finalPath = '%s/maps/bioMart/%s' % ( (settings.getDataDir() if where is None else where), version)
    files = {}
    files['db'] = utils.Acquire(where=where).wget(genQuery()).finalize('%s/data.tsv' % (finalPath))
    return files
  #edef

  def __genQuery(self, url, database, attributes):
    url = self.__versions[version]['url']
    database = self.__versions[version]['database']
    attributes = self.__versions[version]['attributes']

    query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
    query += '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "0" datasetConfigVersion = "0.6" >'
    query += '<Dataset name = "%s" interface = "default" >' % database
    query += ''.join(['<Attribute name = "%s" />' % attr for attr in attributes ])
    query += '</Dataset></Query>'

    return url + query
  #edef

