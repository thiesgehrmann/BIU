import os
import imp

from .. import utils as utils

imp.reload(utils)

class FileManager:

  where    = None
  #fileIndex = None
  str_functions = []
  downloadOnDemand = None

  def __init__(self, where='./', downloadOnDemand=True, **kwargs):
    self.where = os.path.abspath(where)
    self.fileIndex = {}
    #self.fileIndex = self.__urlFileIndex()
    self.str_functions = []
    self.downloadOnDemand = downloadOnDemand
  #edef

  def download(self, what=None, onlyFileNames=False, overwrite=False, **kwargs):
    if what is None:
      what = list(self.fileIndex.keys())
    #fi

    files = {}
    for item in what:
      if item in self.fileIndex:
        url, loc, options = self.fileIndex[item]
        if url is None:
          continue
        #fi
        files[item] = loc
        if not(onlyFileNames):
          utils.downloadIfNotExists(url, loc, overwrite=overwrite, **options)
          if ("tabix" in options) and options["tabix"]:
            utils.tabix(loc, **options)
          #fi
        #fi
      #fi
    #efor
    return files
  #edef

  
  def satisfyRequiredFiles(self, what=None):
    if what is None:
      what = list(self.fileIndex.keys())
    #fi

    for item in what:
      haveItem = self.haveFile(item)
      if haveItem:
        continue
      elif self.downloadOnDemand:
        print("Downloading %s" % item)
        self.download(what=[item])
        if not(self.haveFile(item)):
          print("Error: Failed to download '%s'" % item)
          return False
        #fi
      else:
        return False
      #fi
    #efor
    return True
  #edef

  def haveFile(self, what):
    return os.path.exists(self.getFileName(what))
  #edef

  def getFileName(self, what):
    return self.fileIndex[what][1] if what in self.fileIndex else None
  #edef

  def __urlFileIndex(self):
    files = {}

    return files
  #edef

  def __str__(self):
    dstr  = "%s object\n" % self.__class__.__name__
    dstr += " Where: %s\n" % self.where

    for f in self.str_functions:
      fstr = f(self)
      if fstr is None:
        continue
      #fi
      for line in fstr.split('\n'):
        dstr += ' ' + line + '\n'
    #efor

    dstr += " Files:\n"
    for item in self.fileIndex:
      url, loc, options = self.fileIndex[item]
      dstr += "  * [%s] %s : %s\n" % ('X' if os.path.exists(loc) else ' ', item, loc)
    #efor
    return dstr
  #edef

#eclass
