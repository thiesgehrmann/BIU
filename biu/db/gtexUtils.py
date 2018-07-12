from ..structures import Dataset
from .. import formats
from ..config import settings as settings
from .. import utils

import csv

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

###############################################################################

class GTeX(Dataset):

  __sAttrFieldNames  = [ "SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH", "SMTSPAX", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", "SMAFRZE", "SMGTC", "SME2MPRT", "SMCHMPRS", "SMNTRART", "SMNUMGPS", "SMMAPRT", "SMEXNCRT", "SM550NRM", "SMGNSDTC", "SMUNMPRT", "SM350NRM", "SMRDLGTH", "SMMNCPB", "SME1MMRT", "SMSFLGTH", "SMESTLBS", "SMMPPD", "SMNTERRT", "SMRRNANM", "SMRDTTL", "SMVQCFL", "SMMNCV", "SMTRSCPT", "SMMPPDPR", "SMCGLGTH", "SMGAPPCT", "SMUNPDRD", "SMNTRNRT", "SMMPUNRT", "SMEXPEFF", "SMMPPDUN", "SME2MMRT", "SME2ANTI", "SMALTALG", "SME2SNSE", "SMMFLGTH", "SME1ANTI", "SMSPLTRD", "SMBSMMRT", "SME1SNSE", "SME1PCTS", "SMRRNART", "SME1MPRT", "SMNUM5CD", "SMDPMPRT", "SME2PCTS" ]
  __sPhenoFieldNames = [ "SUBJID", "SEX", "AGE", "DTHHRDY" ]

  versions = {
    "v7" : {
      "gTPM" : 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz',
      "tTPM" : 'GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz',
      "sAttr" : 'GTEx_v7_Annotations_SampleAttributesDS.txt',
      "sPheno" : 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt'
    }
  }

  def __init__(self, version=list(versions.keys())[0], where=None, **kwargs):
    fileIndex = self.__genFileIndex(version, where)
    fm.FileManager.__init__(self, fileIndex, where=where, **kwargs)
    self.version = version
    self._addStrFunction(lambda s: "Note: You must provide your own local copies using the 'localCopy' option!")
    self._addStrFunction(lambda s: "Version: %s" % self.version)

    self._registerObject("_sAttr", pd.read_csv, ["s_attr"], fileIndex["g_tpm"].path, fieldNames=self.__sAttrFieldNames, skiprows=1, delimiter='\t')
    self._registerObject("_sPheno", pd.read_csv, ["s_pheno"], fileIndex["s_pheno"].path, fieldNames=self.__sPhenoFieldNames, skiprows=1, delimiter='\t')
  #edef

  def __genFileIndex(self, version, where=None):
    files = {}
    if where is None:
      where = settings.getDataDir()
    #fi

    files["g_tpm"]   = utils.Acquire('%s/%s' % (where, versions[version]["gTPM"])) 
    files["t_tpm"]   = utils.Acquire('%s/%s' % (where, versions[version]["tTPM"]))
    files["s_attr"]  = utils.Acquire('%s/%s' % (where, versions[version]["sAttr"])) 
    files["s_pheno"] = utils.Acquire('%s/%s' % (where, versions[version]["sPheno"])) 

    return files
  #edef

  ###############################################################################

  def getPersonIDs(self):
    return np.unique(self._sAttr['SAMPID'].apply(lambda x: x.split('-')[1]).values)
  #edef

  def getPersonIDAttrRows(self, personID):
    return self._sAttr[self._sAttr['SAMPID'].apply(lambda x: x.split('-')[1] == personID)]
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
    return np.unique(self._sAttr['SAMPID'])

  def getTissueTypes(self):
    return np.unique(self._sAttr['SMTSD'].values)
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
