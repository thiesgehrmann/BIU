from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

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

  def queryVCF(self, chromosome, start, end):
    return self.query(chromosome, start, end)
  #edef

  def query(self, chromosome, start, end):
    return self.vcf.query(chromosome, start, end)
  #edef

  def queryRegions(self, regions):
    return self.vcf.queryRegions(regions)
  #edef

  def queryCov(self, chromosome, start, end, **kwargs):
    chromosome = str(chromosome)
    return self.cov[chromosome].query(chromosome, start, end, **kwargs)
  #edef

  def queryCovRegions(self, regions):
    R = []
    for (c,s,e) in regions:
      for r in self.queryCov(c, s, e):
        R.append(r)
      #efor
    #efor
    return R
  #edef 
    

#eclass
