
from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

import numpy as np

###############################################################################

versions = { "GRCh37" : {
  "caddURL"      : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz",
  "caddTabixURL" : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi"
  }
}

def urlFileIndex(version):
  files = {}

  files["tsv"] = (versions[version]["caddURL"], 'cadd.tsv.bgz', {})
  files["tsv_tbi"] = (versions[version]["caddTabixURL"], 'cadd.tsv.bgz.tbi', {})

  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class CADD(fm.FileManager):

  version = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["scores"], **kwargs)
    self.version = version

    self.scores = rm.TabixTSVResourceManager(self, "tsv", "tsv_tbi", fieldNames=self._caddFields, **kwargs)

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  #############################################################################

  _caddFields = [ "chrom", "pos", "ref", "alt", "rawscore", "phred" ]

  def query(self, chromosome, start, end=None, alt=None):
    qres = self.scores.query(chromosome, start, (start if (end is None) else end), namedtuple=True)
    if (alt is None) and (end is None):
      resPhred = {}
      for res in qres:
        resPhred[res.alt] = float(res.phred)
      #efor
      return resPhred
    #fi
    if (alt is None) and (end is not None):
      resPhred = {}
      for res in qres:
        resPhred[(int(res.pos),res.alt)] = float(res.phred)
      #efor
      return resPhred
    #fi
    if (alt is not None) and (end is not None):
      resPhred = {}
      for res in [r for r in qres if r.alt == alt]:
        resPhred[int(res.pos)] = float(res.phred)
      #efor
      return resPhred
    #fi
    if (end is None) and (alt is not None):
      relRes = [ r for r in qres if r.alt == alt ]
      if len(relRes) != 1:
        return None
      else:
        return float(relRes[0].phred)
      #fi
    #fi
  #edef

  def queryRegions(self, regions):
    resPhred = {}
    for (c,s,e) in regions:
      qres = self.scores.query(c, s, e, namedtuple=True) 
      for res in qres:
        resPhred[(int(res.pos),res.alt)] = float(res.phred)
      #efor
    #efor
    return resPhred
  #edef

  def regionThresh(self, chromosome, start, end, percentile=95):
    return self.regionsThresh( [(chromosome, start, end)], percentile )
  #edef

  def regionsThresh(self, regions, percentile=95):
    return np.percentile(np.array(list(self.queryRegions(regions).values())), percentile)
  #edef
 

#eclass
   
