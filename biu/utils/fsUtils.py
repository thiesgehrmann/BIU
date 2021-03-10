## File system utilities

import os
import gzip

import glob
import datetime

from . import exeUtils as exe
from . import msgUtils as msg

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

def isEmpty(fileName):
  if not os.path.exists(fileName):
    msg.warning("isEmpty: The file '%s' does not exist." % fileName)
    return True
  elif os.path.isfile(fileName):
    return os.stat(fileName).st_size == 0
  elif os.path.isdir(fileName):
    return len(os.listdir(fileName)) == 0
  else:
    return False
  #fi
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

def most_recent_file(pattern, suffix):
    def datekey(fname):
        return datetime.datetime.strptime(fname.split('.')[-2], '%Y%m%d')
    #edef
    return sorted(glob.glob('%s.*.%s'% (pattern, suffix)), key=datekey)[-1]
#edef