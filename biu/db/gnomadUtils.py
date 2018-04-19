from .. import formats
from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings
from .. import formats

import pandas as pd

###############################################################################

versions = { "GRCh37" : {
  "vcfURL" : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz",
  "vcfTabixURL"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi",
  "covURLProto"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr%s.coverage.txt.gz",
  "chr"           : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" ],
  "covSeqField"   : 1,
  "covBeginField" : 2,
  "covEndField"    : 2
  }
}

def urlFileIndex(version):
  files = {}

  files["vcf"] = (versions[version]["vcfURL"], 'gnomad.vcf.bgz', {})
  files["vcf_tbi"] = (versions[version]["vcfTabixURL"], 'gnomad.vcf.bgz.tbi', {})
  for chrID in versions[version]["chr"]:
    files["chr_%s_cov" % chrID] = (versions[version]["covURLProto"] % chrID, 'gnomad.coverage.chr.%s.tsv.bgz' % chrID, {"tabix": True,
        "bgzip": True,
        "urlIsGzipped":True,
        "seqField" : versions[version]["covSeqField"],
        "beginField" : versions[version]["covBeginField"],
        "endField" : versions[version]["covEndField"]} )
    files["chr_%s_cov_tbi" % chrID] = (None, 'gnomad.coverage.chr.%s.tsv.bgz.tbi' % chrID, {})
  #efor

  return { k : (u, 'gnomad_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class Gnomad(fm.FileManager):

  version = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["vcf"] + [ ("cov", chrID) for chrID in versions[version]["chr" ] ], **kwargs)
    self.version = version

    self.vcf = rm.VCFResourceManager(self, "vcf", "vcf_tbi")
    self.cov = { chrID : rm.TabixTSVResourceManager(self, "chr_%s_cov" % chrID, "chr_%s_cov_tbi" % chrID, fieldNames=self._covEntryFields)
                 for chrID in versions[self.version]["chr" ] }

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  #############################################################################

  _covEntryFields = [ "chrom", "pos", "mean", "median", "q1", "q5", "q10", "q15", "q20", "q25", "q30", "q50", "q100" ]

  def queryVCF(self, *args, **kwargs):
    return self.query(*args, **kwargs)
  #edef

  def query(self, *args, **kwargs):
    return self.queryRegions([ tuple(args) ], **kwargs)
  #edef

  def queryRegions(self, *args, extract=None, sub=None, **kwargs):
    res = list(self.vcf.queryRegions(*args, extract=None, **kwargs))

    #see http://gnomad.broadinstitute.org/faq for the different subpopulations
    if extract  == "summary":
      return self.summary(res, sub=sub)
    #fi
    return formats.VCF(res)
  #edef

  def summary(self, arr, altPos=None, sub=None):
    def allSummary(var, altp):
      gcIndexes  = formats.VCF.genotypeInfoFieldIndexes(altp+1)
      gcmale   = [ var.INFO["GC_Male"][i] for i in gcIndexes ] if "GC_Male" in var.INFO else [0, 0, 0]
      gcfemale   = [ var.INFO["GC_Female"][i] for i in gcIndexes ] if "GC_Female" in var.INFO else [0, 0, 0]

      rr = 0
      r  = 0
      ra = 0
      a  = 0
      aa = 0
      u  = 0

      chrom = var.CHROM.lower()      
      if chrom == 'y':
        a = var.INFO["AC"][altp]
        r = var.INFO["AN"] - a
      elif chrom == 'x':
        rr = gcfemale[0]
        ra = gcfemale[1]
        aa = gcfemale[2]
        r = gcmale[0]
        a = gcmale[1]
      else:
        rr = gcmale[0] + gcfemale[0]
        ra = gcmale[1] + gcfemale[1]
        aa = gcmale[2] + gcfemale[2]
      #fi

      return pd.DataFrame( [(formats.VCF.makeIdentifier(var, altp), rr, r, ra, a, aa, u)],
                           columns = ["id", "RR", "R", "RA", "A", "AA", "O"])
    #edef

    def subSummary(var, altp, sub):
    # See http://gnomad.broadinstitute.org/faq for the different subpopulations
      if sub is None:
        return allSummary(var, altp)
      #fi

      gc = "GC_%s" % sub
      gcIndexes = formats.VCF.genotypeInfoFieldIndexes(altp+1)
      gc        = [ var.INFO[gc][i] for i in gcIndexes ] if gc in var.INFO else [0, 0, 0]
      rr = gc[0]
      r  = 0
      ra = gc[1]
      a  = 0
      aa = gc[2]
      u  = 0
      return pd.DataFrame( [(formats.VCF.makeIdentifier(var, altp), rr, r, ra, a, aa, u)],
                           columns = ["id", "RR", "R", "RA", "A", "AA", "O"])
    #edef

    S = [ subSummary(v, 0 if altPos is None else altPos[i], sub) for i,v in enumerate(arr) ]
    S = pd.concat([ s for s in S if s is not None ])
    return S
  #edef

  def queryCov(self, chromosome, start, end, **kwargs):
    chromosome = str(chromosome)
    return self.cov[chromosome].query(chromosome, start, end, **kwargs)
  #edef

  def queryCovRegions(self, regions, **kwargs):
    R = []
    for (c,s,e) in regions:
      for r in self.queryCov(c, s, e, **kwargs):
        R.append(r)
      #efor
    #efor
    return R
  #edef 
    

#eclass
