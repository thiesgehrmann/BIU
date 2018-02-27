from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

###############################################################################

versions = {
  "human" : "http://geneontology.org/gene-associations/goa_human.gaf.gz",
  "mouse" : "http://geneontology.org/gene-associations/gene_association.mgi.gz",
  "drosophilia" : "http://geneontology.org/gene-associations/gene_association.fb.gz" }

def urlFileIndex(version):
  files = {}

  files["gaf"] = ( versions[version], "geneannots.gaf.gz", {})

  return { k : (u, 'geneOntology_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for version in versions:
    print(" * %s" % version)
  #efor
#edef

###############################################################################

class GO(fm.FileManager):

  _annots = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=[ "_annots" ], **kwargs)
    self.version = version

    self._annots = rm.GAFResourceManager(self, "gaf", skiprows=1, delimiter='\t')

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  ###############################################################################

  def getAnnots(self, objectID):
    return self._annots.getAnnots(objectID)
  #edef

  def getAnnotated(self, annotID):
    return self._annots.getAnnotated(annotID)
  #edef

#eclass
