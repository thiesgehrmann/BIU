import os
import imp
from pathlib import Path

from .. import utils as utils

imp.reload(utils)

class FileManager(object):

  _where    = None
  _fileIndex = None
  _str_functions = []
  _downloadOnDemand = None
  _localCopy = None

  #############################################################################

  def __init__(self, where, fileIndex, objects=None, downloadOnDemand=True, localCopy=None, **kwargs):
    self._where = os.path.abspath(where)
    self._fileIndex = fileIndex
    self._localCopy = localCopy
    self._str_functions = []
    self._downloadOnDemand = downloadOnDemand
    self._makeSymlinks()
    self._objects = objects
  #edef

  #############################################################################

  def _download(self, what=None, overwrite=False, **kwargs):
    if what is None:
      what = list(self._fileIndex.keys())
    #fi

    for item in what:
      if item in self._fileIndex:
        url, fname, options = self._fileIndex[item]
        loc = self.getFileName(item)
        if url is None:
          continue
        #fi

        utils.downloadIfNotExists(url, loc, overwrite=overwrite, **options)
        if ("tabix" in options) and options["tabix"]:
          utils.tabix(loc, **options)
        #fi
      #fi
    #efor
  #edef

  #############################################################################
  
  def satisfyRequiredFiles(self, what=None):
    if what is None:
      what = list(self._fileIndex.keys())
    #fi

    for item in what:
      haveItem = self.haveFile(item)
      if haveItem:
        continue
      elif self._downloadOnDemand:
        print("Downloading %s" % item)
        self._download(what=[item])
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

  #############################################################################

  def haveFile(self, what):
    return os.path.exists(self.getFileName(what))
  #edef

  def getFileName(self, what):
    return '%s/%s' % (self._where, self._fileIndex[what][1]) if what in self._fileIndex else None
  #edef

  #############################################################################

  def _makeSymlinks(self):
    if self._localCopy is not None:
      for item in [ i for i in self._localCopy if (i in self._fileIndex)]:
        loc = self._localCopy[item]
        currLoc = self.getFileName(item)
        if os.path.exists(loc):
          dirname = os.path.dirname(currLoc)
          if not(os.path.exists(dirname)):
            utils.mkdirp(dirname)
          #fi
          if os.path.exists(currLoc) and os.path.islink(currLoc) and (str(Path(currLoc).resolve()) != loc):
            utils.runCommand("unlink '%s'" % currLoc)
            p = utils.runCommand("ln -s '%s' '%s'" % (loc, currLoc))
            print(( "Made symbolic link for '%s'" if p == 0 else "Error using local copy of '%s'") % item)
          elif os.path.exists(currLoc) and os.path.islink(currLoc) and (str(Path(currLoc).resolve()) == loc):
            print("Same symbolic link already exists for '%s'" %item)
          elif not(os.path.exists(currLoc)):
            p = utils.runCommand("ln -s '%s' '%s'" % (loc, currLoc))
            print(( "Made symbolic link for '%s'" if p == 0 else "Error using local copy of '%s'") % item)
          else:
            print("Could not use local copy of '%s' as file already exists at '%s'" % (item, currLoc))
          #fi
        else:
          print("Error: Local copy of '%s' does not exist" % item)
        #fi
      #efor
    #fi
  #edef

  #############################################################################

  def addStrFunction(self, func):
    self._str_functions.append(func)
  #edef

  def __str__(self):
    dstr  = "%s object\n" % self.__class__.__name__
    dstr += " Where: %s\n" % self._where

    for f in self._str_functions:
      fstr = f(self)
      if fstr is None:
        continue
      #fi
      for line in fstr.split('\n'):
        dstr += ' ' + line + '\n'
    #efor

    if self._objects is not None:
      dstr += " Objects:\n"
      for o in self._objects:
        if isinstance(o, tuple):
          oname, idx = o
          dstr += "  * [%s] %s[%s]\n" % ('X' if getattr(self,oname)[idx].__initialized__() else ' ', oname, str(idx))
        else:
          dstr += "  * [%s] %s\n" % ('X' if getattr(self, o).__initialized__() else ' ', o)
        #fi
      #efor
    #fi

    dstr += " Files:\n"
    for item in self._fileIndex:
      loc = self.getFileName(item)
      if os.path.islink(loc):
        dstr += "  * [%s] %s : %s -> %s\n" % ('S' if os.path.exists(loc) else ' ', item, loc, Path(loc).resolve())
      else:
        dstr += "  * [%s] %s : %s\n" % ('X' if os.path.exists(loc) else ' ', item, loc)
      #fi
    #efor
    return dstr
  #edef

  #############################################################################

#eclass
