from ..structures import Dataset
from .. import formats
from ..config import settings as settings
from .. import utils as utils

pd = utils.py.loadExternalModule("pandas")

import os

###############################################################################


class HAGR(Dataset):

  version = None
  versions = {
    "current": {
      "humanGenes"    : ("http://genomics.senescence.info/genes/human_genes.zip", 'genage_human.csv', ','),
      "modelGenes"    : ("http://genomics.senescence.info/genes/models_genes.zip", "genage_models.csv", ','),
      "longevityMap"  : ("http://genomics.senescence.info/longevity/longevity_genes.zip", "longevity.csv", ','),
      "drugage"       : ("http://genomics.senescence.info/drugs/dataset.zip", "drugage.csv", ','),
      "anage"         : ("http://genomics.senescence.info/species/dataset.zip", "anage_data.txt", '\t') }
  }

  def __init__(self, version=list(versions.keys())[0], where=None):
    fileIndex = self.__genFileIndex(version, where=where)
    Dataset.__init__(self, fileIndex)
    self.version = version
    for oname in fileIndex:
      self._registerObject(oname, pd.read_csv, [ oname ], fileIndex[oname].path, delimiter=self.versions[version][oname][2])
    #efor

  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/HAGR_%s' % ( (settings.getWhere() if where is None else where), version)
     files = {}
     for oname in self.versions[version]:
       url, selFileName, delim = self.versions[version][oname]
       files[oname] = utils.Acquire(where=where).curl(url).unzip(selFileName).finalize('%s/%s' % (finalPath, selFileName))
     #efor
     return files
  #edef

#eclass
