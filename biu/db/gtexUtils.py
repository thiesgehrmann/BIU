import csv
import os
from . import fileManager as fm
from .. import utils
import imp
import pandas as pd

imp.reload(fm)
imp.reload(utils)
from collections import namedtuple

###############################################################################

versions = {
  "v7" : {
    "gTPM" : None,
    "tTPM" : None,
    "sAttr" : None,
    "sPheno" : None
  }

}

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
    fm.FileManager.__init__(self, where, **kwargs)
    self.version = version
    self.fileIndex = self.__urlFileIndex()
    self.str_functions.append(lambda s: "Note: You must provide your own local copies using the 'localCopy' option!")
    self.str_functions.append(lambda s: "Version: %s" % self.version)

    def loadedObjects(s):
      dstr = 'Objects:\n'
      dstr += " * [%s] _gTPMSource\n" % ('X' if s._gTPMSource is not None else ' ')
      dstr += " * [%s] _tTPMSource\n" % ('X' if s._tTPMSource is not None else ' ')
      dstr += " * [%s] _sAttrSource\n" % ('X' if s._sAttrSource is not None else ' ')
      dstr += " * [%s] _sPhenoSource\n" % ('X' if s._sPhenoSource is not None else ' ')
      return dstr
    #edef
    self.str_functions.append(lambda s: loadedObjects(self))
  #edef

  def __urlFileIndex(self):
    files = {}

    files["g_tpm"] = (versions[self.version]["gTPM"], self.where + '/genes_TPM.tsv.gz', {})
    files["t_tpm"] = (versions[self.version]["tTPM"], self.where + '/transcripts_TPM.tsv.gz', {})
    files["s_attr"] = (versions[self.version]["gTPM"], self.where + '/sample_attributes.tsv', {})
    files["s_pheno"] = (versions[self.version]["gTPM"], self.where + '/sample_phenotypes.tsv', {})
    return files
  #edef

  _gTPMSource = None
  _tTPMSource = None
  _sAttrSource = None
  _sPhenoSource = None

  _sAttrFieldNames = [ "SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN", "SMTS", "SMTSD", "SMUBRID", "SMTSISCH", "SMTSPAX", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", "SMAFRZE", "SMGTC", "SME2MPRT", "SMCHMPRS", "SMNTRART", "SMNUMGPS", "SMMAPRT", "SMEXNCRT", "SM550NRM", "SMGNSDTC", "SMUNMPRT", "SM350NRM", "SMRDLGTH", "SMMNCPB", "SME1MMRT", "SMSFLGTH", "SMESTLBS", "SMMPPD", "SMNTERRT", "SMRRNANM", "SMRDTTL", "SMVQCFL", "SMMNCV", "SMTRSCPT", "SMMPPDPR", "SMCGLGTH", "SMGAPPCT", "SMUNPDRD", "SMNTRNRT", "SMMPUNRT", "SMEXPEFF", "SMMPPDUN", "SME2MMRT", "SME2ANTI", "SMALTALG", "SME2SNSE", "SMMFLGTH", "SME1ANTI", "SMSPLTRD", "SMBSMMRT", "SME1SNSE", "SME1PCTS", "SMRRNART", "SME1MPRT", "SMNUM5CD", "SMDPMPRT", "SME2PCTS" ]
  def _requiresAttrSource(self):
    if self._sAttrSource is None:
      if not(self.satisfyRequiredFiles(["s_attr"])):
        return False
      #fi
      #self._sAttrSource = utils.readNamedColumnTsvFile(self.getFileName("s_attr"), columnNames = self._sAttrFieldNames, delimiter='\t', skip=1)
      self._sAttrSource = pd.read_csv(self.getFileName("s_attr"), delimiter='\t', names=self._sAttrFieldNames, skiprows=1)
    #fi
    return True
  #edef


#eclass
