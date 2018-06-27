# Utilities to acquire files

from . import fsUtils as fs
from . import msgUtils as msg
from . import exeUtils as exe
from ..config import settings

import os
from collections import namedtuple
import hashlib

## Tdo

###############################################################################

acquireStep = namedtuple("AcquireStep", [ "step", "pargs", "kwargs" ])

class Acquire(object):

  __slots__ = [ '__fileName', '__finalName', '__steps', '__redo', '__dlDir' ]

  def __init__(self, fileName=None, finalName=None, steps=[], redo=False, where=None):
    self.__fileName  = fileName
    self.__finalName = fileName if finalName is None else finalName
    self.__steps     = steps
    self.__redo      = redo
    self.__dlDir     = settings.getDownloadDir() if where is None else where
  #edef

  def __str__(self):
    dstr = "Acquire object.\n"
    dstr += ' Re-do steps: %s\n' % ('yes' if self.__redo else 'no')
    dstr += " Current steps:\n"
    for step in self.__steps:
      dstr += '  * %s\n' % step.step
    #efor
    return dstr
  #edef

  #############################################################################

  @property
  def exists(self):
    if self.__finalName is None:
      return False
    #fi
    return (os.path.exists(self.__finalName) and self.__checkExistsTag(self.__finalName))
  #edef

  @property
  def path(self):
    return self.__finalName
  #edef

  #############################################################################

  def acquire(self):
    if self.exists and not self.__redo:
      return self.__finalName
    #fi

    for step in self.__steps:
      oldFileName = self.__fileName
      status = getattr(self, '_' + step.step)(*step.pargs, **step.kwargs)
      print(status)
      if status != 0:
        msg.error("Could not complete step '%s'" % step.step)
        self.__rmExistsTag(self.__fileName)
        return None
      else:
        self.__setExistsTag(self.__fileName) # Set exist tag.
      #fi
    #efor
    self.__finalName = self.__fileName
    return self.__fileName
  #edef

  #############################################################################

  def redo(self, redo=True):
    self.__redo = redo
    return self
  #edef

  #############################################################################

  def __addStep(self, step, finalName=None):
    newSteps = self.__steps + [step]
    if finalName is None:
      finalName = self.__finalName
    #fi
    return Acquire(fileName=self.__fileName, finalName=finalName, steps=newSteps, redo=self.__redo, where=self.__dlDir)
  #edef

  def curl(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("curl", pargs, kwargs))
  def lftp(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("lftp", pargs, kwargs))
  def local(self, fileName, *pargs, **kwargs):
    return self.__addStep(acquireStep("local", tuple([fileName]) + pargs, kwargs), finalName=fileName)
  def wget(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("wget", pargs, kwargs))
  def touch(self, fileName=None, *pargs, **kwargs):
    if fileName is None:
      from datetime import datetime
      fileName = self.__dlDir + '/touchedFile.' + str(self.__downloadHash(str(datetime.now())))
    #fi
    return self.__addStep(acquireStep("touch", tuple([fileName]) + pargs, kwargs), finalName=fileName)
  def cmd(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("cmd", pargs, kwargs))
  def call(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("call", pargs, kwargs))
  def cat(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("cat", pargs, kwargs))
  def ls(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("ls", pargs, kwargs))
  def unzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("unzip", pargs, kwargs))
  def bunzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("bunzip", pargs, kwargs))
  def gunzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("gunzip", pargs, kwargs))
  def untar(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("untar", pargs, kwargs))
  def select(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("select", pargs, kwargs))
  def sort(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("sort", pargs, kwargs))
  def tabix(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("tabix", pargs, kwargs))
  def gzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("gzip", pargs, kwargs))
  def bgzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("bgzip", pargs, kwargs))
  def bzip(self, *pargs, **kwargs):
    return self.__addStep(acquireStep("bzip", pargs, kwargs))
  def finalize(self, finalName, *pargs, **kwargs):
    return self.__addStep(acquireStep("finalize", tuple([finalName]) + pargs, kwargs), finalName=finalName)

  #############################################################################  

  def __downloadHash(self, arguments):
    sha1 = hashlib.sha1()
    for a in arguments:
      if a is not None:
        sha1.update(a.encode('utf-8'))
      #fi
    #efor
    return sha1.hexdigest()
  #edef

  def __getExtension(self, fileName):
    return ''
  #edef

  def __checkExistsTag(self, fileName):
    return os.path.exists(fileName + '.__exists__')
  #edef

  def __setExistsTag(self, fileName):
    fs.touchFile(fileName + '.__exists__')
  #edef

  def __rmExistsTag(self, fileName):
    fs.rmFile(fileName + '.__exists__')
  #edef

  def _curl(self, url, cookieURL=None, username=None, password=None, ext=None):
    ext = self.__getExtension(url) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ url, cookieURL, username, password ])
    self.__fileName = self.__dlDir + '/' + curlHash +  ext
    print(self.__fileName)

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    fs.mkdirname(self.__fileName)
    cookieFile = None
    if (cookieURL is not None) or (username is not None) or (password is not None):
      cookieFile = self.__fileName + '.cookie'
      p = exe.runCommand('curl -b %s --user "%s:%s" "%s"' % (cookieFile, username, password, cookieURL))
      if p != 0:
        msg.error("Could not set cookie... for site '%s'" % cookieURL)
        return p
      #fi
    #fi
    p = exe.runCommand("curl -L %s '%s' > '%s'" % ( '-c "%s"' % cookieFile if cookieFile is not None else  '', url, self.__fileName), shell=True, verbose=True)
    return p
  #edef

  def _lftp(self, server, location, username, password, ext=None):
    ext = self.__getExtension(location) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ server, location, username, password ])
    self.__fileName = self.__dlDir + '/' + curlHash
    print(self.__fileName)

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    fs.mkdirname(self.__fileName)
    cmd = "echo -en  'open \"%s\"\\nuser \"%s\" \"%s\"\\ncat \"%s\"' | lftp > '%s'" % (server, username, password, location, self.__fileName)
    p = exe.runCommand(cmd, shell=True, verbose=True)
    return p
  #edefi

  def _wget(self, url, ext=None):
    ext = self.__getExtension(url) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ url ])
    self.__fileName = self.__dlDir + '/' + curlHash +  ext

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    cmd = "wget -O '%s' '%s'" % ( self.__fileName, url )
    p = exe.runCommand(cmd, verbose=True)
    return p
  #edef

  def _local(self, fileName):
    self.__fileName = fileName
    if os.path.isfile(self.__fileName):
      return 0
    else:
      return 1
    #fi
  #edef

  def _touch(self, fileName):
    self.__fileName = fileName
    return fs.touchFile(fileName)
  #edef

  def _call(self, cmd):
    if '%s' in cmd:
      cmd = cmd % (self.__fileName)
    else:
      cmd = "%s '%s'" % (cmd, self.__fileName)
    #fi
    print(exe.getCommandOutput(cmd, verbose=True, shell=True).decode('utf-8'))
    return 0
  #edef

  def _cmd(self, cmd):
    oldFile = self.__fileName
    self.__fileName = self.__fileName + '.cmd'

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi
    if '%s' in cmd:
      cmd = cmd % oldFile
    else:
      cmd = cmd + ' ' + oldFile
    #fi
    p = exe.runCommand("%s > '%s'" % (cmd, self.__fileName), shell=True, verbose=True)
    return p
  #edef

  def _cat(self):
    print(exe.getCommandOutput("cat '%s'" % self.__fileName).decode('utf-8'))
    return 0
  #edef

  def _ls(self):
    print(str(exe.getCommandOutput("ls -R '%s'" % self.__fileName).decode('utf-8')))
    return 0
  #edef

  def _unzip(self, fileName=None):
    zipFileName = self.__fileName
    outDirName  = '%s.unzipped' % zipFileName
    self.__fileName = outDirName
    if fileName is not None:
      self.__fileName = self.__fileName + '/' + fileName
    #fi

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("unzip -o -d '%s' '%s'" % (outDirName, zipFileName), verbose=True)
    return p
  #edef

  def _bunzip(self):
    return 0
  #edef

  def _gunzip(self):
    gzipFileName = self.__fileName
    outDirName  = '%s.gunzipped' % gzipFileName
    self.__fileName = outDirName

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("gunzip < '%s' > '%s'" % (gzipFileName, self.__fileName), verbose=True, shell=True)
    return p
  #edef

  def _untar(self, fileName=None):
    tarFile = self.__fileName
    outDirName = tarFile + '.untar'
    self.__fileName = outDirName
    if fileName is not None:
      self.__fileName = outDirName + '/' + fileName
    #fi

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("tar -xf '%s' -C '%s'" % (tarFile, outDirName), verbose=True)
    return p
  #edef

  def _select(self, fileName):
    self.__fileName = self.__fileName + '/' + fileName
  #edef

  def _tabix(self, fileType=None, seq=0, start=1, end=2):
    bgzipFileName = self.__fileName
    outName = bgzipFileName + '.tbi'
    self.__fileName = outName

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    if fileType is not None:
      cmd = "tabix -p %s '%s'" % (fileType, bgzipFileName)
    else:
      cmd = "tabix -s %d -b %d -e %d %s" % (seq, start, end, bgzipFileName)
    #fi
    p = exe.runCommand(cmd, verbose=True)
    return p
  #edef

  def _sort(self, options=None):
    fileName = self.__fileName
    outName = fileName + '.sorted'
    self.__fileName = outName

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("sort %s < '%s' > '%s'" % ('' if options is None else options, fileName, outName), shell=True)
    return p
  #edef

  def _bgzip(self):
    oldFile = self.__fileName
    self.__fileName = oldFile + '.bgz'

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("bgzip < '%s' > '%s'" % (oldFile, self.__fileName), verbose=True, shell=True)
    return p
  #edef

  def _gzip(self):
    oldFile = self.__fileName
    self.__fileName = oldFile + '.gz'

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("gzip < '%s' > '%s'" % (oldFile, self.__fileName), verbose=True, shell=True)
    return p
  #edef

  def _bzip(self):
    oldFile = self.__fileName
    self.__fileName = oldFile + '.bz'

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    p = exe.runCommand("bzip < '%s' > '%s'" % (oldFile, self.__fileName), verbose=True, shell=True)
    return p
  #edef

  def _finalize(self, fileName, ln=False):
    oldFile = self.__fileName
    self.__fileName = os.path.realpath(os.path.expanduser(fileName))

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    fs.mkdirname(self.__fileName)

    if ln:
      p = exe.runCommand("ln -s '%s' '%s'" % (oldFile, self.__fileName), verbose=True)
    else:
      p = exe.runCommand("cp '%s' '%s'" % (oldFile, self.__fileName), verbose=True)
    #fi
    return p
  #edef

#eclass
 
