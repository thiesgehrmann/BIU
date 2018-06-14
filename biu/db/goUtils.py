from ..structures import Dataset
from .. import formats
from .. import utils
from ..config import settings

###############################################################################

class GO(Dataset):

  versions = {
    "human" : "http://geneontology.org/gene-associations/goa_human.gaf.gz",
    "mouse" : "http://geneontology.org/gene-associations/gene_association.mgi.gz",
    "drosophilia" : "http://geneontology.org/gene-associations/gene_association.fb.gz" }

  def __init__(self, version=list(versions.keys())[0], where=None, **kwargs):
    fileIndex = self.__genFileIndex(version, where=where)
    Dataset.__init__(self, fileIndex, **kwargs)
    self.version = version

    self._registerObject('gaf', formats.GAF, ['gaf'], fileIndex["gaf"].path, skiprows=1, delimiter='\t')

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/geneOntology_%s' % ( (settings.getWhere() if where is None else where), version)
     url = self.versions[version]
     files = {}
     files['gaf'] = utils.Acquire(where=where).curl(url).gunzip().finalize('%s/annots.gaf' % finalPath)
     return files
  #edef


  ###############################################################################

  @property
  def gaf(self):
    return self._getObject("gaf")
  #edef

  @property
  def annotations(self):
    return self._getObject("gaf").annotations
  #edef

  @property
  def objects(self):
    return self._getObject("gaf").objects
  #edef

  def getAnnots(self, objectID):
    return self._getObject("gaf").getAnnots(objectID)
  #edef

  def getAnnotated(self, annotID):
    return self._getObject("gaf").getAnnotated(annotID)
  #edef

  def enrich(self, *pargs, **kwargs):
    return self._getObject("gaf").enrich(*pargs, **kwargs)
  #edef

#eclass
