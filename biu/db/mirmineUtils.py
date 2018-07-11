from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings
from .. import utils as utils

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

import os

## NOTE NOTE NOTE
## NOTE: Here there is a dirty hack because we have ZIP files as a source, which our model is currently unable to handle...


###############################################################################

versions = { "current" : [] }

def urlFileIndex():

  files = {}

  files["miRmine.zip"]   = ("http://guanlab.ccmb.med.umich.edu/mirmine/miRmine.zip", "miRmine.zip", {})

  files["cellLines"]   = (None, "miRmine.zip.unpacked/miRmine-cell-lines.xlsx", {})
  files["tissues"]   = (None, "miRmine.zip.unpacked/miRmine-tissues.xlsx", {})
  files["info"]   = (None, "miRmine.zip.unpacked/miRmine-info.txt", {})

  return { k : (u, 'miRmine/%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef


class MiRmine(fm.FileManager):

  version = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=["_tissues", "_cellLines", "_info"], **kwargs)
    if not(os.path.exists('%s.unpacked' % self.getFileName("miRmine.zip"))):
      self._download([ "miRmine.zip" ])
      utils.runCommand("unzip -d '%s' '%s'" % ( '%s.unpacked' % self.getFileName("miRmine.zip"), self.getFileName("miRmine.zip")))
    #fi
   
    self._cellLines = rm.XLSXResourceManager(self, "cellLines")
    self._tissues  = rm.XLSXResourceManager(self, "tissues")
    self._info     = rm.TSVResourceManager(self, "info")
  #edef

  def _loadCellLinesDF(self):
    if not(self._cellLines.isSheetLoaded("Sheet1")):
      self._cellLines.loadSheet("Sheet1", header=True)
    #fi
  #edef

  def _loadTissuesDF(self):
    if not(self._tissues.isSheetLoaded("Sheet1")):
      self._tissues.loadSheet("Sheet1", header=True)
    #fi
  #edef

  def _getAccessionType(self, accession):
    res = self._info[self._info["Experiment Accession"] == accession].values

    if len(res) != 1:
      return None
    else:
      tissue, cellLine = res[0][1:3]
      if pd.isnull(tissue):
        return "C"
      else:
        return "T"
      #fi
    #fi
  #edef

  def _getAccessionColName(self, accession):
    res = self._info[self._info["Experiment Accession"] == accession].values

    if len(res) != 1:
      return None
    else:
      tissue, cellLine = res[0][1:3]
      if pd.isnull(tissue):
        return "%s (%s)" % (accession, cellLine)
      else:
        return "%s (%s)" % (accession, tissue)
      #fi
    #fi
  #edef

  def getAccessions(self):
    return list(set(self._info["Experiment Accession"].values))
  #edef

  def getExpr(self, accessions):

    if isinstance(accessions, str):
      accessions = [ accessions ]
    #fi

    matureID, precursorID = self._getmiRNANameColumns()
    D = {}
    D["matureID"] = matureID
    D["precursorID"] = precursorID
    for a in accessions:
      D[a] = self._getTissueExpr(a) if self._getAccessionType(a) == 'T' else self._getCellLineExpr(a)
    #efor
    return pd.DataFrame(data=D)[ [ "matureID", "precursorID"] + accessions ]
  #edef

  def _getmiRNANameColumns(self):
    if self._cellLines.isSheetLoaded("Sheet1"):
      return list(zip(*self._cellLines["Sheet1"][["Mature miRNA ID", "Precursor miRNA ID"]].values))
    elif self._tissues.isSheetLoaded("Sheet1"):
      return list(zip(*self._tissues["Sheet1"][["Mature miRNA ID", "Precursor miRNA ID"]].values))
    else:
      self._loadCellLinesDF()
      return list(zip(*self._cellLines["Sheet1"][["Mature miRNA ID", "Precursor miRNA ID"]].values))
    #fi
  #edef

  def _getCellLineExpr(self, accession):
    self._loadCellLinesDF()
    aColName = self._getAccessionColName(accession)
    return self._cellLines["Sheet1"][aColName].values
  #edef

  def _getTissueExpr(self, accession):
    self._loadTissuesDF()
    aColName = self._getAccessionColName(accession)
    return self._tissues["Sheet1"][aColName].values
  #edef

#eclass
