from ..structures import Dataset
from ..config import settings as settings
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule("pandas")

###############################################################################

versions = { "GRCh37" : {
  "vcfURL" : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz",
  "vcfTabixURL"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi",
  "covURLProto"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr%s.coverage.txt.gz",
  "chr"           : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" ],
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

class Gnomad(Dataset):

  versions = { "GRCh37" : {
    "vcf"           : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz",
    "tbi"           : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi",
    "covURLProto"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr%s.coverage.txt.gz",
    "chr"           : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" ],
    "covSeqField"   : 1,
    "covBeginField" : 2,
    "covEndField"    : 2
    }
  }

  version = None

  def __init__(self, version=list(versions.keys())[0], where=None):
    fileIndex = self.__genFileIndex(version, where)
    Dataset.__init__(self, fileIndex)
    self.version = version

    self._registerObject("_vcf", formats.VCF, [ "vcf", "tbi" ], fileIndex["vcf"].path, tabix=True)
    covEntryFields = [ "chrom", "pos", "mean", "median", "q1", "q5", "q10", "q15", "q20", "q25", "q30", "q50", "q100" ]
    for chrID in self.versions[self.version]["chr"]:
      self._registerObject("_cov_%s" % chrID, formats.Tabix, [ "cov_%s" % chrID, "cov_%s_tbi" % chrID ], fileIndex["cov_%s" % chrID].path, fieldNames = covEntryFields)
    #efor

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/gnomad_%s' % ( (settings.getWhere() if where is None else where), version)
     vData = self.versions[version]
     files = {}
     files['vcf'] = utils.Acquire(where=where).curl(vData["vcf"]).finalize('%s/gnomad.vcf.bgz' % finalPath)
     files['tbi'] = utils.Acquire(where=where).curl(vData["tbi"]).finalize('%s/gnomad.vcf.bgz.tbi' % finalPath)
     for chrID in vData["chr"]:
       covTXT = utils.Acquire(where=where).curl(vData["covURLProto"] % chrID).gunzip().bgzip()
       files['cov_%s' % chrID] = covTXT.finalize('%s/gnomad.coverage.chr.%s.tsv.bgz' % (finalPath, chrID))
       files['cov_%s_tbi' % chrID] = covTXT.tabix(seq=vData["covSeqField"], begin=vData["covBeginField"], end=vData["covEndField"]).finalize('%s/gnomad.coverage.chr.%s.tsv.bgz.tbi' % (finalPath, chrID))
     return files
  #edef

  #############################################################################

  def queryVCF(self, *args, **kwargs):
    return self.query(*args, **kwargs)
  #edef

  def query(self, *args, **kwargs):
    return self.queryRegions([ tuple(args) ], **kwargs)
  #edef

  def queryRegions(self, *args, extract=None, sub=None, **kwargs):
    res = self._vcf.queryRegions(*args, extract="raw", **kwargs)

    #see http://gnomad.broadinstitute.org/faq for the different subpopulations
    if extract  == "summary":
      return self.summary(res, sub=sub)
    #fi
    return formats.VCF(res)
  #edef

  def summary(self, arr, altPos=None, sub=None):
    def allSummary(var, altp):
      gcIndexes  = formats.VCF.genotypeInfoFieldIndexes(altp)
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
        a = var.INFO["AC"][altp-1]
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

      return pd.DataFrame( [(formats.VCF.makeIdentifier(var, altp-1), rr, r, ra, a, aa, u)],
                           columns = ["id", "RR", "R", "RA", "A", "AA", "O"])
    #edef

    def subSummary(var, altp, sub):
    # See http://gnomad.broadinstitute.org/faq for the different subpopulations
      if sub is None:
        return allSummary(var, altp)
      #fi

      gc = "GC_%s" % sub
      gcIndexes = formats.VCF.genotypeInfoFieldIndexes(altp)
      gc        = [ var.INFO[gc][i] for i in gcIndexes ] if gc in var.INFO else [0, 0, 0]
      rr = gc[0]
      r  = 0
      ra = gc[1]
      a  = 0
      aa = gc[2]
      u  = 0
      return pd.DataFrame( [(formats.VCF.makeIdentifier(var, altp-1), rr, r, ra, a, aa, u)],
                           columns = ["id", "RR", "R", "RA", "A", "AA", "O"])
    #edef

    S = [ subSummary(v, 0 if altPos is None else altPos[i], sub) for i,v in enumerate(arr) ]
    S = pd.concat([ s for s in S if s is not None ])
    return S
  #edef

  def queryCov(self, chromosome, start, end, **kwargs):
    chromosome = str(chromosome)
    oname = "_cov_%s" % chromosome
    return self._getObject(oname).query(chromosome, start, end, **kwargs)
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
    

  def getVar(self, *pargs, **kwargs):
    return self._vcf.getVar(*pargs, **kwargs)
  #edef

  def whoHas(self, *pargs, **kwargs):
    return self._vcf.whoHas(*pargs, **kwargs)
  #edef

#eclass
