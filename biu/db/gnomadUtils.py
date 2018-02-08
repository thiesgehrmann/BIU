from collections import namedtuple
import tabix

import vcf as vcf

from . import fileManager as fm
from .. import utils as utils
import imp
imp.reload(fm)
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

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class Gnomad(fm.FileManager):

  version = None

  def __init__(self, version=list(versions.keys())[0], where='./', **kwargs):
    fm.FileManager.__init__(self, where, **kwargs)
    self.version = version
    self.fileIndex = self.__urlFileIndex()
    self.str_functions.append(lambda s: "Version: %s" % self.version)

    def loadedObjects(s):
      dstr = 'Objects:\n'
      dstr += " * [%s] _vcfSource\n" % ('X' if s._vcfSource is not None else ' ')
      dstr += " * [%s] _covSource\n" % ('X' if s._covSource is not None else ' ')
      return dstr
    #edef

    self.str_functions.append(lambda s: loadedObjects(s))
  #edef

  def __urlFileIndex(self):
    files = {}

    files["vcf"] = (versions[self.version]["vcfURL"], self.where + '/gnomad.vcf.bgz', {})
    files["vcf_tbi"] = (versions[self.version]["vcfTabixURL"], self.where + '/gnomad.vcf.bgz.tbi', {})
    for chrID in versions[self.version]["chr"]:
      files["chr_%s_cov" % chrID] = (versions[self.version]["covURLProto"] % chrID, self.where + '/gnomad.coverage.chr.%s.tsv.bgz' % chrID, {"tabix": True, 
          "bgzip": True,
          "urlIsGzipped":True,
          "seqField" : versions[self.version]["covSeqField"],
          "beginField" : versions[self.version]["covBeginField"],
          "endField" : versions[self.version]["covEndField"]} )
      files["chr_%s_cov_tbi" % chrID] = (None, self.where + '/gnomad.coverage.chr.%s.tsv.bgz.tbi' % chrID, {})
    #efor

    return files
  #edef

  #############################################################################

  _vcfSource = None
  _covSource = None

  def _requireVcfSource(self):
    if self._vcfSource is None:
      if not(self.satisfyRequiredFiles(["vcf", "vcf_tbi"])):
        print("Could not satisfy files")
        return False
      #fi
      self._vcfSource = vcf.Reader(filename=self.getFileName("vcf"), compressed=True)
    #fi
    return True
  #edef

  def _requireCovSource(self, chromosome):
    if self._covSource is None:
      self._covSource = { chrID: None for chrID in versions[self.version]["chr"] }
    #fi
    if self._covSource[chromosome] is None:
      if not(self.satisfyRequiredFiles(["chr_%s_cov" % chromosome, "chr_%s_cov_tbi" % chromosome])):
        return False
      #fi
      self._covSource[chromosome] = tabix.open(self.getFileName("chr_%s_cov" % chromosome))
    #fi
    return True
  #edef

  def queryVCF(self, chromosome, start, end):
    if not(self._requireVcfSource()):
      print("Could not initiate vcf source")
      return None
    #fi
    return utils.vcfQueryWrapper(self._vcfSource, chromosome, start, end)
  #edef

  covEntryFields = [ "chrom", "pos", "mean", "median", "q1", "q5", "q10", "q15", "q20", "q25", "q30", "q50", "q100" ]
  covEntry = namedtuple("gnomadCoverageEntry", covEntryFields)

  def queryCov(self, chromosome, start, end):
    chromosome = str(chromosome)
    if not(self._requireCovSource(chromosome)):
      return None
    #fi
    return [ self.covEntry(*r) for r in utils.tabixQueryWrapper(self._covSource[chromosome], chromosome, start, end) ]
  #edef
    

#eclass

class gnomadLookup:

  _vcfSource = None
  _covSource = None

  def __init__(self, _vcfSource, _covSource):
    self._vcfSource = vcf.Reader(file=_vcfSource)
    self._covSource = tabix.open(_covSource)
  #edef

  def queryVCF(self, chromosome, start, end):
    self._vcfSource.fetch(chromosome, start, end)
    return res
  #edef

  def query(self, chromosome, start, end):
    return self.queryVCF(chromosome, start, end)
  #edef

  def queryCov(self, chromosome, start, end):
    return self._covSource.query(chromosome, start, end)
  #edef

#eclass
    
