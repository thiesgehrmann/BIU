from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings
from .. import utils as utils

import os

## NOTE NOTE NOTE
## NOTE: Here there is a dirty hack because we have ZIP files as a source, which our model is currently unable to handle...


###############################################################################

versions = { "current" : [] }

def urlFileIndex():

  files = {}

  files["human_genes.zip"]   = ("http://genomics.senescence.info/genes/human_genes.zip", 'human_genes.zip', {})
  files["model_genes.zip"]   = ("http://genomics.senescence.info/genes/models_genes.zip", "model_genes.zip", {})
  files["longevity_map.zip"] = ("http://genomics.senescence.info/longevity/longevity_genes.zip", "longevity_map.zip", {})
  files["drugage.zip"]       = ("http://genomics.senescence.info/drugs/dataset.zip", "drugage.zip", {})
  files["anage.zip"]         = ("http://genomics.senescence.info/species/dataset.zip", "anage.zip", {})

  files["human_genes"]   = (None, "human_genes.zip.unpacked/genage_human.csv", {})
  files["model_genes"]   = (None, "model_genes.zip.unpacked/genage_models.csv", {})
  files["longevity_map"] = (None, "longevity_map.zip.unpacked/longevity.csv", {})
  files["drugage"]       = (None, "drugage.zip.unpacked/drugage.csv", {})
  files["anage"]         = (None, "anage.zip.unpacked/anage_data.txt", {})

  return { k : (u, 'hagr_%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class HAGR(fm.FileManager):

  version = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=["human_genes", "model_genes", "longevity_genes", "drug_age", "an_age"], **kwargs)
    for what in [ "human_genes.zip", "model_genes.zip", "longevity_map.zip", "drugage.zip", "anage.zip" ]:
      self._download([ what ])
      if not(os.path.exists('%s.unpacked' % self.getFileName(what))):
        utils.runCommand("unzip -d '%s' '%s'" % ( '%s.unpacked' % self.getFileName(what), self.getFileName(what)))
      #fi
    self.human_genes     = rm.TSVResourceManager(self, "human_genes", delimiter=',', **kwargs)
    self.model_genes     = rm.TSVResourceManager(self, "model_genes", delimiter=',', **kwargs)
    self.longevity_genes = rm.TSVResourceManager(self, "longevity_map", delimiter=',', **kwargs)
    self.drug_age        = rm.TSVResourceManager(self, "drugage", delimiter=',', **kwargs)
    self.an_age          = rm.TSVResourceManager(self, "anage", delimiter='\t', **kwargs)
  #edef

#eclass
