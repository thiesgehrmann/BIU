from ..structures import Dataset
from .. import formats
from ..config import settings
from .. import utils

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

class ClinVar(Dataset):

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

  def __init__(self, version=list(versions.keys())[0], where=None):
    fileIndex = self.__genFileIndex(version)
    Dataset.__init__(self, fileIndex)
    self.version = version

    self._registerObject("_vcf", formats.VCF, [ "vcf", "vcf_tbi" ], fileIndex["vcf"].path, tabix=True)
    self._registerObject("_summary", formats.Tabix, ["sum", "sum_tbi" ], fileIndex["sum"].path, fieldNames=clinVarSummaryFields)

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/clinvar_%s' % ( (settings.getWhere() if where is None else where), version)
     vData = self.versions[version]
     files = {}
     files['vcf']     = utils.Acquire(where=where).curl(vData["vcfURL"]).finalize('%s/clinVar.vcf.bgz' % finalPath)
     files['vcf_tbi'] = utils.Acquire(where=where).curl(vData["vcfTbiURL"]).finalize('%s/clinVar.vcf.bgz.tbi' % finalPath)
     summaryAcq = utils.Acquire().curl(vData["summaryURL"])\
                                 .gunzip()\
                                 .sort("-t $'\\t' -k19,19V -k 20,21n")\
                                 .cmd("awk -F $'\\t' 'BEGIN {OFS = FS} { if($19 != \"na\"){ print $0}}'")\
                                 .bgzip()
     files['sum']     = summaryAcq.finalize("%s/summary.tsv.bgz" % finalPath)
     files['sum_tbi'] = summaryAcq.tabix(seq=19, start=20, end=21).finalize("%s/summary.tsv.bgz.tbi" % finalPath)

     return files
  #edef

  ###############################################################################

  def queryVCF(self, *args, **kwargs):
    return self.query(*args, **kwargs)
  #edef

  def queryRegions(self, *args, **kwargs):
    return self._vcf.queryRegions(*args, **kwargs)
  #edef

  def query(self, chromosome, start, end, **kwargs):
    return self.queryRegions([(chromosome, start, end)], **kwargs)
  #edef

  def querySummary(self, chromosome, start, end, namedtuple=True, **kwargs):
    return self._summary.query(chromosome, start, end, namedtuple=namedtuple, **kwargs)
  #edef

#eclass
###############################################################################
