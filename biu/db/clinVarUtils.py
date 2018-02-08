import csv
import os
import tabix
from . import fileManager as fm
from .. import utils
import imp

imp.reload(fm)
imp.reload(utils)
from collections import namedtuple
import vcf

###############################################################################

versions = {
  "GRCh37" :
  { "vcfURL"     : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    "vcfTbiURL"  : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi",
    "summaryURL" : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
    "seqField"   : 19,
    "beginField" : 20,
    "endField"   : 21,
    "delim"      : '\t' },
  "GRCh38":
  { "vcfURL"     : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
    "vcfTbiURL"  : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi",
    "summaryURL" : "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
    "seqField"   : 19,
    "beginField" : 20,
    "endField"   : 21,
    "delim"      : '\t' }
}

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

clinVarSummaryFields = [ "alleleid",
                         "type",
                         "name",
                         "geneid",
                         "genesymbol",
                         "hgnc_id",
                         "clinicalsignificance",
                         "clinsigsimple",
                         "lastevaluated",
                         "rs",
                         "nsv_esv",
                         "rcvaccession",
                         "phenotypeids",
                         "phenotypelist",
                         "origin",
                         "originsimple",
                         "assembly",
                         "chromosomeaccession",
                         "chromosome",
                         "start",
                         "stop",
                         "referenceallele",
                         "alternateallele",
                         "cytogenetic",
                         "reviewstatus",
                         "numbersubmitters",
                         "guidelines",
                         "testedingtr",
                         "otherids",
                         "submittercategories" ]

CVS = namedtuple("clinVarSummary", clinVarSummaryFields);

###############################################################################

def listVersions():
  print("Available versions:")
  for version in versions:
    print(" * %s" % version)
  #efor
#edef

###############################################################################

class ClinVar(fm.FileManager):

  genomeID = None
  where    = None
  fileIndex = None

  def __init__(self, version=list(versions.keys())[0], where='./', **kwargs):
    fm.FileManager.__init__(self, where, **kwargs)
    self.version = version
    self.fileIndex = self.__urlFileIndex()
    self.str_functions.append(lambda s: "Version: %s" % self.version)

    def loadedObjects(s):
      dstr = 'Objects:\n'
      dstr += " * [%s] _vcfSource\n" % ('X' if s._vcfSource is not None else ' ')
      dstr += " * [%s] _summarySource\n" % ('X' if s._summarySource is not None else ' ')
      return dstr
    #edef
    self.str_functions.append(lambda s: loadedObjects(self))
  #edef

  def __urlFileIndex(self):
    files = {}

    files["vcf"] = (versions[self.version]["vcfURL"], self.where + '/clinVar.vcf.bgz', {})
    files["vcf_tbi"] = (versions[self.version]["vcfTbiURL"], self.where + '/clinVar.vcf.bgz.tbi', {})
    files["sum"] = (versions[self.version]["summaryURL"], self.where + '/summary.tsv.bgz', {"bgzip":True,
        "tabix":True,
        "seqField":19,
        "beginField":20,
        "endField":21,
        "inject":"sort -t $'\\t' -k19,19V -k 20,21n | awk -F $'\\t' 'BEGIN {OFS = FS} { if($19 != \"na\"){ print $0}}'"})
    files["sum_tbi"] = (None, self.where + '/summary.tsv.bgz.tbi', {})
    
    return files
  #edef

  ###############################################################################

  _vcfSource = None
  _summarySource = None

  def _requireVcfSource(self):
    if self._vcfSource is None:
      if not(self.satisfyRequiredFiles(["vcf", "vcf_tbi"])):
        return False
      #fi
      self._vcfSource = vcf.Reader(filename=self.getFileName("vcf"), compressed=True)
    #fi
    return True
  #edef
  def _requireSummarySource(self):
    if self._summarySource is None:
      if not(self.satisfyRequiredFiles(["sum", "sum_tbi"])):
        return False
      #fi
      self._summarySource = tabix.open(self.getFileName("sum"))
    #fi
    return True
  #fi

  def queryVCF(self, chromosome, start, end):
    if not(self._requireVcfSource()):
      return None
    #fi
    return utils.vcfQueryWrapper(self._vcfSource, chromosome, start, end)
  #edef

  def querySummary(self, chromosome, start, end):
    if not(self._requireSummarySource()):
      return None
    #fi
    res = [ CVS(*r) for r in utils.tabixQueryWrapper(self._summarySource, chromosome, start, end) ]
    return res
  #edef

#eclass

#eclass
###############################################################################

def castClinVarSummaryRow(row):
  return row;
#edef

###############################################################################

clinVarDelim = '\t'

def readClinVarSummary(fileName, setTypes=False, delimiter=clinVarDelim):

  cvs = []

  with utils.gzopen(fileName, "r") as ifd:
    for row in csv.reader(ifd, delimiter=delimiter):
      if (len(row) == 0) or (len(row[0]) == 0):
        continue
      elif row[0][0] == "#":
        continue
      if len(row) == len(clinVarSummaryFields):
        if setTypes:
          row = castClinVarSummaryRow(row)
        #fi
        cvs.append(CVS(*row))
      #fi
    #efor
  #ewith

  return cvs

#edef

###############################################################################

def writeClinVarSummary(cvs, fileName, delimiter=clinVarDelim):
  isGzipped = fileName[-2:] == "gz"

  with utils.gzopen(fileName, "w") as ofd:
    ofd.write('#' + delimiter.join(clinVarSummaryFields) + '\n')
    for summary in cvs:
      ofd.write(delimiter.join([ str(x) for x in summary ]) + '\n')
    #efor
  #ewith
#edef

###############################################################################
