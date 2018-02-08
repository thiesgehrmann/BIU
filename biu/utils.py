import gzip
import pathlib
import urllib
import os
import shlex, subprocess, os;

###############################################################################

def mkdirp(directory):
  pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
#edef

###############################################################################

def download(url, fileName, urlIsGzipped=None, fileIsGzipped=None, bgzip=None, inject=None, **kwargs):
  urlIsGzipped  = (url[-2:] == "gz") if urlIsGzipped is None else urlIsGzipped
  fileIsGzipped = (fileName[-2:] == "gz") if fileIsGzipped is None else fileIsGzipped

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
  fullCommand = "curl -L '%s'%s%s%s > '%s'" % (url, 
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
    download(url, fileName, urlIsGzipped, fileIsGzipped, bgzip=bgzip, **kwargs)
    return True
  #fi
  return False
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

def tabixQueryWrapper(struct, chrom, start, end):
  try:
    return struct.query(str(chrom), int(start), int(end))
  except Exception as e:
    print(e)
    return []
  #etry
#edef

###############################################################################

def vcfQueryWrapper(struct, chrom, start, end):
  try:
    return struct.fetch(str(chrom), int(start), int(end))
  except Exception as e:
    print(e)
    return []
  #etry
#edef

###############################################################################

def gzopen(fileName, mode):
  isGzipped = fileName[-2:] == "gz"
  readMode  = mode == "r"

  if isGzipped and readMode:
    return gzip.open(fileName, "rb")
  elif isGzipped and not(readMode):
    return gzip.open(fileName, "wb")
  elif not(isGzipped) and readMode:
    return open(fileName, "r")
  else:
    return open(fileName, "w")
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
