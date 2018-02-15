
from . import fileManager as fm
from . import resourceManager as rm

###############################################################################

versions = {
  "v7" : {
    "gTPM" : None,
    "tTPM" : None,
    "sAttr" : None,
    "sPheno" : None
  }

}

def urlFileIndex(version):
  files = {}
  files["g_tpm"] = (versions[version]["gTPM"], 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz', {})
  files["t_tpm"] = (versions[version]["tTPM"], 'GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz', {})
  files["s_attr"] = (versions[version]["s_attr"], 'GTEx_v7_Annotations_SampleAttributesDS.txt', {})
  files["s_pheno"] = (versions[version]["s_pheno"], 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt', {})
  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

class GTeX(fm.FileManager):

  genomeID = None
  where    = None
  fileIndex = None

  def __init__(self, version=list(versions.keys())[0], where='./', **kwargs):
    fm.FileManager.__init__(self, where, urlFileIndexversion(), ["sAttr", "sPheno", "gTPM", "tTPM"], **kwargs)
    self.version = version
    self.addStrFunction(lambda s: "Note: You must provide your own local copies using the 'localCopy' option!")
    self.addStrFunction(lambda s: "Version: %s" % self.version)

    self.sAttr  = rm.TSVResourceManager(self, "s_attr", fieldNames = self._sAttrFieldNames, skiprows=1)
    self.sPheno = rm.TSVResourceManager(self, "s_pheno", fieldNames = self._sPhenoFieldNames, skiprows=1)
    self.gTPM   = rm.TSVResourceManager(self, "g_tpm", skiprows=0)
    self.tTPM   = rm.TSVResourceManager(self, "t_tpm", skiprows=0)
  #edef

  _sAttrFieldNames  = [ "SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH", "SMTSPAX", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", "SMAFRZE", "SMGTC", "SME2MPRT", "SMCHMPRS", "SMNTRART", "SMNUMGPS", "SMMAPRT", "SMEXNCRT", "SM550NRM", "SMGNSDTC", "SMUNMPRT", "SM350NRM", "SMRDLGTH", "SMMNCPB", "SME1MMRT", "SMSFLGTH", "SMESTLBS", "SMMPPD", "SMNTERRT", "SMRRNANM", "SMRDTTL", "SMVQCFL", "SMMNCV", "SMTRSCPT", "SMMPPDPR", "SMCGLGTH", "SMGAPPCT", "SMUNPDRD", "SMNTRNRT", "SMMPUNRT", "SMEXPEFF", "SMMPPDUN", "SME2MMRT", "SME2ANTI", "SMALTALG", "SME2SNSE", "SMMFLGTH", "SME1ANTI", "SMSPLTRD", "SMBSMMRT", "SME1SNSE", "SME1PCTS", "SMRRNART", "SME1MPRT", "SMNUM5CD", "SMDPMPRT", "SME2PCTS" ]
  _sPhenoFieldNames = [ "SUBJID", "SEX", "AGE", "DTHHRDY" ]

#eclass
