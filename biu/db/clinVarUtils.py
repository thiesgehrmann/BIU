from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

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


###############################################################################

def urlFileIndex(version):
  files = {}

  files["vcf"]     = (versions[version]["vcfURL"], 'clinVar.vcf.bgz', {})
  files["vcf_tbi"] = (versions[version]["vcfTbiURL"], 'clinVar.vcf.bgz.tbi', {})
  files["sum"]     = (versions[version]["summaryURL"], 'summary.tsv.bgz', {"bgzip":True,
      "tabix":True,
      "seqField":19,
      "beginField":20,
      "endField":21,
      "inject":"sort -t $'\\t' -k19,19V -k 20,21n | awk -F $'\\t' 'BEGIN {OFS = FS} { if($19 != \"na\"){ print $0}}'"})
  files["sum_tbi"] = (None, 'summary.tsv.bgz.tbi', {})

  return { k : (u, 'clinvar_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for version in versions:
    print(" * %s" % version)
  #efor
#edef

###############################################################################

class ClinVar(fm.FileManager):

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=[ "summary", "vcf" ], **kwargs)
    self.version = version

    # Define the objects in the clinVar fields
    self.summary = rm.TabixTSVResourceManager(self, "sum", "sum_tbi", fieldNames = clinVarSummaryFields)
    self.vcf     = rm.VCFResourceManager(self, "vcf", "vcf_tbi")

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  ###############################################################################

  def queryVCF(self, chromosome, start, end):
    return self.vcf.query(chromosome, start, end)
  #edef

  def querySummary(self, chromosome, start, end):
    return self.summary.query(chromosome, start, end, namedtuple=True)
  #edef

#eclass

#eclass
###############################################################################
