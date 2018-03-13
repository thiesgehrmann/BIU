from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

import itertools

###############################################################################

versions = { "current":
  { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "M", "X", "Y" ] }
}

def urlFileIndex(version):
  files = {}

  chrs = versions[version]["chrs"]
  for chrID in chrs:
    files["vcf_%s" % chrID] = (None, "tbx/merged.chr%s.vcf.bgz" % chrID, {})
    files["vcf_%s_tbi" % chrID] = (None, "tbx/merged.chr%s.vcf.bgz.tbi" % chrID, {})
  #efor

  files["phen"] = (None, "phen218.txt", {})

  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class LLS(fm.FileManager):

  version = None
  vcf = None

  def __init__(self, version=list(versions.keys())[0], where="/exports/molepi/LLSSEQ", **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=[ ("vcf", chrID) for chrID in versions[version]["chrs"] ], where=where, **kwargs)
    self.version = version

    self.vcf = {}
    for chrID in versions[self.version]["chrs"]:
      self.vcf[chrID] = rm.VCFResourceManager(self, "vcf_%s" % chrID, "vcf_%s_tbi" % chrID)
    #efor

    self.phenotypes = rm.TSVResourceManager(self, "phen", delimiter=' ') 

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def query(self, chrID, start, end):
    chrID = str(chrID)
    if chrID in self.vcf:
      return self.vcf[chrID].query(chrID, start, end)
    else:
      utils.error("Could not find chromosome '%s'" % chrID)
      return iter(())
    #fi
  #edef

  def queryRegions(self, regions):
    R = []
    for (c,s,e) in regions:
      for r in self.query(c, s, e):
        R.append(r)
      #efor
    #efor
    return R
  #edef

#eclass
