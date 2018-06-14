## File system utilities

import os
import gzip

from . import exeUtils as exe

###############################################################################

def mkdirname(fileName):
  dirname = os.path.dirname(fileName)
  return mkdirp(dirname)
#edef

###############################################################################

def mkdirp(directory):
  #pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
  return exe.runCommand("mkdir -p '%s'" % directory)
#edef

###############################################################################

def rmFile(fileName):
  if os.path.isfile(fileName):
    return exe.runCommand("rm '%s'" % fileName)
  #fi
  return 1
#edef

###############################################################################

def touchFile(fileName):
  mkdirname(fileName)
  p = exe.runCommand("touch '%s'" % fileName)
  return p
#edef

###############################################################################

def gzopen(fileName, mode="r", gzmode='t', **kwargs):
  isGzipped = fileName[-2:] == "gz"
  if isGzipped:
    return gzip.open(fileName, mode + gzmode, **kwargs)
  else:
    return open(fileName, mode)
  #fi
#edef

###############################################################################
