from ..structures import Dataset
from ..config import settings as settings
from .. import formats
from .. import utils

import itertools

###############################################################################

def urlFileIndex(version):
  files = {}

  chrs = versions[version]["chrs"]
  for chrID in chrs:
    files["vcf_%s" % chrID] = (None, "tbx/merged.bbmri.chr%s.vcf.bgz" % chrID, {})
    files["vcf_%s_tbi" % chrID] = (None, "tbx/merged.bbmri.chr%s.vcf.bgz.tbi" % chrID, {})
  #efor

  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class BBMRI(Dataset):

  versions = { "current":
    { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "M", "X", "Y" ] }
  }

  version = None
  vcf = None

  def __init__(self, version=list(versions.keys())[0], where="/exports/molepi/BBMRISEQ", **kwargs):
    fileIndex = self.__genFileIndex(version, where)
    Dataset.__init__(self, fileIndex)
    self.version = version

    
    for chrID in self.versions[self.version]["chrs"]:
      self._registerObject('vcf_%s' % chrID, formats.VCF, [ "vcf_%s" % chrID, "vcf_%s_tbi" % chrID ], fileIndex["vcf_%s" % chrID].path, tabix=True)
    #efor

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     files = {}
     for chrID in self.versions[version]["chrs"]:
       files['vcf_%s' % chrID] = utils.Acquire("%s/tbx/merged.bbmri.chr%s.vcf.bgz" % (where, chrID), where=where)
       files['vcf_%s_tbi' % chrID] = utils.Acquire("%s/tbx/merged.bbmri.chr%s.vcf.bgz.tbi" % (where, chrID), where=where)
     #efor

     return files
  #edef

  #############################################################################

  def query(self, chrID, start, end, **kwargs):
    chrID = str(chrID)
    oname = "vcf_%s" % chrID
    if self._objectExists(oname):
      return self._getObject(oname).query(chrID, start, end, **kwargs)
    else:
      utils.error("Could not find chromosome '%s'" % chrID)
      return iter(())
    #fi
  #edef

  def queryRegions(self, regions, extract=None, **kwargs):
    R = []
    for (c,s,e) in regions:
      R.extend(list(self.query(c,s,e, **kwargs)))
    #efor
    return formats.VCF.extract(R, extract=extract)
  #edef

  def getVar(self, chromosome, *pargs, **kwargs):
    chromosome = str(chromosome)
    oname = "vcf_%s" % chromosome
    if not self._objectExists(oname):
      utils.error("Could not find chromosome '%s'" % chromosome)
      return None
    #fi
    return self._getObject(oname).getVar(chromosome, *pargs, **kwargs)
  #edef

  def whoHas(self, chromosome, *pargs, **kwargs):
    chromosome = str(chromosome)
    oname = "vcf_%s" % chromosome
    if not self._objectExists(oname):
      utils.error("Could not find chromosome '%s'" % chromosome)
      return None
    #fi
    return self._getObject(oname).whoHas(chromosome, *pargs, **kwargs)
  #edef

#eclass
