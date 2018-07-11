
from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

from .. import utils

import csv

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

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
  # We don't add anything to the filenames for GTeX, because it is behind this ugly google login.
  # We don't want people to have to rename the files or figure out the correct filenames
  files = {}
  files["g_tpm"] = (versions[version]["gTPM"], 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz', {})
  files["t_tpm"] = (versions[version]["tTPM"], 'GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz', {})
  files["s_attr"] = (versions[version]["sAttr"], 'GTEx_v7_Annotations_SampleAttributesDS.txt', {})
  files["s_pheno"] = (versions[version]["sPheno"], 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt', {})
  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class GTeX(fm.FileManager):

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["sAttr", "sPheno"], **kwargs)
    self.version = version
    self.addStrFunction(lambda s: "Note: You must provide your own local copies using the 'localCopy' option!")
    self.addStrFunction(lambda s: "Version: %s" % self.version)

    self.sAttr  = rm.TSVResourceManager(self, "s_attr", fieldNames = self._sAttrFieldNames, skiprows=1)
    self.sPheno = rm.TSVResourceManager(self, "s_pheno", fieldNames = self._sPhenoFieldNames, skiprows=1)

  #edef

  _sAttrFieldNames  = [ "SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH", "SMTSPAX", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", "SMAFRZE", "SMGTC", "SME2MPRT", "SMCHMPRS", "SMNTRART", "SMNUMGPS", "SMMAPRT", "SMEXNCRT", "SM550NRM", "SMGNSDTC", "SMUNMPRT", "SM350NRM", "SMRDLGTH", "SMMNCPB", "SME1MMRT", "SMSFLGTH", "SMESTLBS", "SMMPPD", "SMNTERRT", "SMRRNANM", "SMRDTTL", "SMVQCFL", "SMMNCV", "SMTRSCPT", "SMMPPDPR", "SMCGLGTH", "SMGAPPCT", "SMUNPDRD", "SMNTRNRT", "SMMPUNRT", "SMEXPEFF", "SMMPPDUN", "SME2MMRT", "SME2ANTI", "SMALTALG", "SME2SNSE", "SMMFLGTH", "SME1ANTI", "SMSPLTRD", "SMBSMMRT", "SME1SNSE", "SME1PCTS", "SMRRNART", "SME1MPRT", "SMNUM5CD", "SMDPMPRT", "SME2PCTS" ]
  _sPhenoFieldNames = [ "SUBJID", "SEX", "AGE", "DTHHRDY" ]

  ###############################################################################

  def getPersonIDs(self):
    return np.unique(self.sAttr['SAMPID'].apply(lambda x: x.split('-')[1]).values)
  #edef

  def getPersonIDAttrRows(self, personID):
    return self.sAttr[self.sAttr['SAMPID'].apply(lambda x: x.split('-')[1] == personID)]
  #edef

  def getPersonIDSamples(self, personID, smafrze="RNASEQ"):
    personRows = self.getPersonIDAttrRows(personID)
    return personRows[personRows['SMAFRZE'] == smafrze]['SAMPID'].values
  #edef

  def getPersonIDTissueSampleID(self, personID, tissueType):
    personRows = self.getPersonIDAttrRows(personID)
    return personRows[personRows['SMTSD'] == tissueType]['SAMPID'].values
  #edef

  def getSampleIDs(self):
    return np.unique(self.sAttr['SAMPID'])

  def getTissueTypes(self):
    return np.unique(self.sAttr['SMTSD'].values)
  #edef

  ###############################################################################
  # Query the gene expression level database

  def getGeneExpr(self, sampleIDs):
    if not(self.satisfyRequiredFiles(["g_tpm"])):
      utils.error("Could not satisfy dependency")
      return None
    #fi

    if len(set(sampleIDs) & set(self.getSampleIDs())) != len(set(sampleIDs)):
      utils.error("The sampleIDs you specified are not present in GTEx")
      return None
    #fi

    return pd.read_csv(self.getFileName("g_tpm"), usecols=["Name"] + list(sampleIDs), delimiter='\t', skiprows=2)
  #edef

  def getTranscriptExpr(self, sampleIDs):
    if not(self.satisfyRequiredFiles(["t_tpm"])):
      utils.error("Could not satisfy dependency")
      return None
    #fi

    if len(set(sampleIDs) & set(self.getSampleIDs())) != len(set(sampleIDs)):
      utils.error("The sampleIDs you specified are not present in GTEx")
      return None
    #fi

    return pd.read_csv(self.getFileName("t_tpm"), usecols=["transcript_id"] + list(sampleIDs), delimiter='\t', skiprows=2)
  #edef
#eclass
