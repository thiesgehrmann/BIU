from .. import utils

pd = utils.py.loadExternalModule("pandas")
xl = utils.py.loadExternalModule('xlrd')

class XLS(object):

  _fileName = None
  _xlBook = None
  _data = None

  def __init__(self, fileName):
    self._fileName = fileName
    self._xlBook = xl.open_workbook(fileName)
    self._data = {}
  #edef

  def loadSheet(self, sheetName, index=None, columns=None, header=0, names=None, skiprows=0, refresh=False):
    if sheetName not in self._xlBook.sheet_names():
      utils.error("Sheet '%s' not in workbook '%s'." % (sheetName, self._fileName))
      return None
    #fi

    sheetData = pd.read_excel(io=self._xlBook, sheet_name=sheetName, header=header, skiprows=skiprows, index=index, engine="xlrd")

    self._data[sheetName] = sheetData
    return self._data[sheetName]
  #edef

  def isSheetLoaded(self, sheetName):
    return sheetName in self._data
  #edef

  def __getitem__(self, sheetName):
    if (sheetName not in self._data):
      return self.loadSheet(sheetName)
    else:
      return self._data[sheetName]
    #fi
  #edef

  def __str__(self):
    dstr  = "XLS object\n"
    dstr += " Where: %s\n" % (self._fileName if self._fileName is not None else hex(id(self)))
    dstr += " Sheets:\n"
    for sheet in self._xlBook.sheet_names():
      dstr += "  * [%s] %s\n" % ( 'X' if sheet in self._data else ' ', sheet)
    #efor
    return dstr
  #edef 

  def toTSV(self, sheetName, outFile, **kwargs):
    self[sheetName].to_csv(outFile, **kwargs)
  #edef

#eclass
