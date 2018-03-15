import gzip
import pathlib
import urllib
import os
import shlex, subprocess, os;
import csv
import sys
import hashlib
from collections import namedtuple

from ..config import settings as settings

###############################################################################

def stripkwargs(kwargs):
  internalKwargs = [ "where", "version", "objects" ]

  return { k : kwargs[k] for k in kwargs if k not in internalKwargs }
#edef

###############################################################################

def mkdirname(fileName):
  dirname = os.path.dirname(fileName)
  return mkdirp(dirname)
#edef

###############################################################################

def mkdirp(directory):
  #pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
  return runCommand("mkdir -p '%s'" % directory)
#edef

###############################################################################

def rmFile(fileName):
  if os.path.exists(fileName):
    if os.path.isfile(fileName):
      return runCommand("rm '%s'" % fileName)
    #fi
  #fi
#edef

###############################################################################

def touchFile(fileName):
  mkdirname(fileName)
  p = runCommand("touch '%s'" % fileName)
  return p
#edef

###############################################################################

def dbm(message):
  if settings.getDebugState():
    for line in str(message).split('\n'):
      if settings.getDebugStream == 'stdout':
        sys.stdout.write('D: %s\n' % line)
      else:
        sys.stderr.write('D: %s\n' % line)
      #fi
    #efor
  #fi
#edef

###############################################################################

def error(message):
  for line in str(message).split('\n'):
    sys.stderr.write('E: %s\n' % line)
  #efor
#edef

###############################################################################

def warning(message):
  for line in str(message).split('\n'):
    sys.stderr.write('W: %s\n' % line)
  #efor
#edef

###############################################################################

def download(url, fileName, urlIsGzipped=None, fileIsGzipped=None, bgzip=None, zipPath=None, inject=None, curlCommand=None, **kwargs):
  """
    Available options
      *urlIsGzipped: Boolean, is the URL gZipped? (Detected from url if not provided)
      *fileIsGzipped: Boolean, should the file be gZipped? (detected from filename if not provided)
      *bgzip: Boolean, should the output file be bgzipped, for use with tabix?
      *inject: String, command to perform after curl, possible unzipping, but before possible gzipping/bgzipping
      *curlCommand: String, command to perform in place of curl, if there is a strange URL.
  """
  urlIsGzipped  = (url[-2:] == "gz") if ((urlIsGzipped is None) and (url is not None)) else urlIsGzipped
  fileIsGzipped = (fileName[-2:] == "gz") if fileIsGzipped is None else fileIsGzipped

  #urlIsZipped   = (url[-3:] == 'zip') if urlIsZipped is None else urlIsZipped
  #fileIzZipped  = (fileName[-3:] == 'zip') if fileIsZipped is None else fileIsZipped 


    # Make the directory if it doesn't exist yer
  mkdirp(os.path.dirname(fileName))

  pre  = None
  post = None
  if urlIsGzipped and bgzip:
    pre = "zcat"
    post = "bgzip"
  elif bgzip:
    post = "bgzip"
  elif (urlIsGzipped == fileIsGzipped):
    pass
  elif (urlIsGzipped): # The URL is gzipped, but the output is NOT
    pre = "zcat"
  else: # The URL is not gzipped, but the output IS
    post = "gzip"
  #fi

  if curlCommand is None:
    curlCommand = "curl -L '%s'" % url
  #fi

  fullCommand = "%s%s%s%s > '%s'" % (curlCommand, 
    (' | %s' % pre) if not(pre is None) else '',
    (' | %s' % inject) if not(inject is None) else '',
    (' | %s' % post) if not(post is None) else '',
    fileName)
  #print(fullCommand)
  p = runCommand(fullCommand, shell=True)

  return p

#edef

###############################################################################

def downloadIfNotExists(url, fileName, urlIsGzipped=None, fileIsGzipped=None, overwrite=False, bgzip=None, **kwargs):
  if not(os.path.isfile(fileName)) or overwrite:
    return download(url, fileName, urlIsGzipped, fileIsGzipped, bgzip=bgzip, **kwargs)
  #fi
  return 0
#edef

###############################################################################

def bgzip(fileName, out):
  runCommand("cat '%s' | bgzip > '%s'" % (fileName, out), shell=True)
#edef

def tabix(fileName, seqField=1, beginField=2, endField=3, ignoreExtension=False, **kwargs):
  if ignoreExtension or (fileName[-3:] != "bgz"):
    error("File '%s' does not appear to be bgzipped. Set ignoreExtension=True if you want to cirumvent this")
    return None
  else:
    #print("tabix -s %d -b %d -e %d %s" % (seqField, beginField, endField, fileName))
    r = runCommand("tabix -s %d -b %d -e %d %s" % (seqField, beginField, endField, fileName))
    return fileName + ".tbi"
  #fi
#edef

###############################################################################

def tabixQueryWrapper(struct, chrom, start, end):
  try:
    return struct.query(str(chrom), int(start), int(end))
  except Exception as e:
    print(e)
    return []
  #etry
#edef

###############################################################################

def gzopen(fileName, mode="r", **kwargs):
  isGzipped = fileName[-2:] == "gz"
  if isGzipped:
    return gzip.open(fileName, mode, **kwargs)
  else:
    return open(fileName, mode, **kwargs)
  #fi
#edef


###############################################################################


def runCommand(cmd, bg=False, stdin=None, stdout=None, stderr=None, shell=False, verbose=False):
  if verbose:
    print(cmd)
  #fi

  if not shell:
    cmd = shlex.split(cmd)
  #fi

  p = subprocess.Popen(cmd, stdin=stdin, stdout=stdout, stderr=stderr, shell=shell);
  if bg:
    return p;
  else:
    (pid, r) = os.waitpid(p.pid, 0);
    return r;
  #fi
#edef

###############################################################################

def getCommandOutput(cmd, stderr=None, shell=False, verbose=False):
  if verbose:
    print(cmd)
  #fi

  if not shell:
    cmd = shlex.split(cmd)
  #fi

  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell);
  (pid, r) = os.waitpid(p.pid, 0);
  return p.communicate()
#edef

###############################################################################

def readNamedColumnTsvFile(fileName,
                           columnNames, 
                           namedTupleName="NamedColumnFile", 
                           delimiter='\t', 
                           quotechar=None, 
                           skip=0, 
                           maxLines=None, 
                           columnTypes=None, 
                           **kwargs):

  rowType = namedtuple(namedTupleName, columnNames)
  nFields = len(columnNames)
  D = []
  i = 0
  with open(fileName, "r") as ifd:
    reader = csv.reader(ifd, delimiter=delimiter, quotechar=quotechar)
    for row in reader:
      i += 1
      if (i < skip):
        continue
      elif (maxLines is not None) and (i > maxLines):
        break
      elif len(row) != nFields:
        print("Error: line %d does not have %d fields!" % (i, nFields))
      else:
        D.append( rowType(*row) )
      #fi
    #efor
  #ewith
  return D
#edef


###############################################################################

def hashArray(arr, strategy="tmb", f=hashlib.md5):
  h = f()

  if strategy == 'tmb':
    middle = int(len(arr) / 2)
    for o in arr[:10] + arr[middle:middle+10] + arr[:-10]:
      h.update(str(o).encode())
    #efor
  #fi

  return h.hexdigest()
#edef

###############################################################################
