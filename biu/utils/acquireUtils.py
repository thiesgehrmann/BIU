# Utilities to acquire files

from . import fsUtils as fs
from . import msgUtils as msg
from . import exeUtils as exe
from ..config import settings

import os
from collections import namedtuple
import hashlib
import ftplib

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
    """
    exists: True if the final file exists. False otherwise
    """
    if self.__finalName is None:
      return False
    #fi
    if len(self.__steps) == 0: # The file specified was local
      return os.path.exists(self.__finalName)
    #fi
    return (os.path.exists(self.__finalName) and self.__checkExistsTag(self.__finalName))
  #edef

  @property
  def path(self):
    """
    path: The path to the final file. Only exists if finalize() exists in the pipeline
    """
    return self.__finalName
  #edef

  @property
  def steps(self):
    """
    steps: The steps in the pipeline
    """
    return self.__steps
  #edef

  #############################################################################

  def acquire(self):
    """
    acquire: Execute the acquire pipeline
    """
    if self.exists and not self.__redo:
      return self.__finalName
    #fi

    fs.mkdirp(self.__dlDir)

    for step in self.__steps:
      oldFileName = self.__fileName
      status = getattr(self, '_' + step.step)(*step.pargs, **step.kwargs)
      if status != 0:
        msg.error("Could not complete step '%s'" % step.step)
        self.__rmExistsTag(self.__fileName)
        raise RuntimeError("Could not acquire this file. Failed at step '%s'." % step.step)
      else:
        self.__setExistsTag(self.__fileName) # Set exist tag.
      #fi
    #efor
    self.__finalName = self.__fileName
    return self.__fileName
  #edef

  #############################################################################

  def redo(self, redo=True):
    """
    redo: Set the redo flag
    Inputs: redo : Boolean, redo the whole pipeline or not. (default True) """
    self.__redo = redo
    return self
  #edef

  def where(self, where):
    """
    where: Set the download directory
    Inputs: where : Path to download directory
    """
    self.__dlDir = where
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
    """ See help for _curl """
    return self.__addStep(acquireStep("curl", pargs, kwargs))
  def ftp(self, *pargs, **kwargs):
    """ See help for _ftp """
    return self.__addStep(acquireStep("ftp", pargs, kwargs))
  def lftp(self, *pargs, **kwargs):
    """ See help for _lftp """
    return self.__addStep(acquireStep("lftp", pargs, kwargs))
  def local(self, fileName, *pargs, **kwargs):
    """ See help for _local """
    return self.__addStep(acquireStep("local", tuple([fileName]) + pargs, kwargs), finalName=fileName)
  def wget(self, *pargs, **kwargs):
    """ See help for _wget """
    return self.__addStep(acquireStep("wget", pargs, kwargs))
  def touch(self, fileName=None, *pargs, **kwargs):
    """ See help for _touch """
    if fileName is None:
      from datetime import datetime
      fileName = self.__dlDir + '/touchedFile.' + str(self.__downloadHash(str(datetime.now())))
    #fi
    return self.__addStep(acquireStep("touch", tuple([fileName]) + pargs, kwargs), finalName=fileName)
  def merge(self, *pargs, **kwargs):
    """ See help for _merge """
    return self.__addStep(acquireStep("merge", pargs, kwargs))
  def cmd(self, *pargs, **kwargs):
    """ See help for _cmd """
    return self.__addStep(acquireStep("cmd", pargs, kwargs))
  def func(self, *pargs, **kwargs):
    """ See help for _func """
    return self.__addStep(acquireStep("func", pargs, kwargs))
  def call(self, *pargs, **kwargs):
    """ See help for _call """
    return self.__addStep(acquireStep("call", pargs, kwargs))
  def cat(self, *pargs, **kwargs):
    """ See help for _cat """
    return self.__addStep(acquireStep("cat", pargs, kwargs))
  def ls(self, *pargs, **kwargs):
    """ See help for _ls """
    return self.__addStep(acquireStep("ls", pargs, kwargs))
  def unzip(self, *pargs, **kwargs):
    """ See help for _unzip """
    return self.__addStep(acquireStep("unzip", pargs, kwargs))
  def bunzip(self, *pargs, **kwargs):
    """ See help for _bunzip """
    return self.__addStep(acquireStep("bunzip", pargs, kwargs))
  def gunzip(self, *pargs, **kwargs):
    """ See help for _gunzip """
    return self.__addStep(acquireStep("gunzip", pargs, kwargs))
  def untar(self, *pargs, **kwargs):
    """ See help for _untar """
    return self.__addStep(acquireStep("untar", pargs, kwargs))
  def select(self, *pargs, **kwargs):
    """ See help for _select """
    return self.__addStep(acquireStep("select", pargs, kwargs))
  def sort(self, *pargs, **kwargs):
    """ See help for _sort """
    return self.__addStep(acquireStep("sort", pargs, kwargs))
  def tabix(self, *pargs, **kwargs):
    """ See help for _tabix """
    return self.__addStep(acquireStep("tabix", pargs, kwargs))
  def gzip(self, *pargs, **kwargs):
    """ See help for _gzip """
    return self.__addStep(acquireStep("gzip", pargs, kwargs))
  def bgzip(self, *pargs, **kwargs):
    """ See help for _bgzip """
    return self.__addStep(acquireStep("bgzip", pargs, kwargs))
  def bzip(self, *pargs, **kwargs):
    """ See help for _bzip """
    return self.__addStep(acquireStep("bzip", pargs, kwargs))
  def finalize(self, finalName, *pargs, **kwargs):
    """ See help for _finalize """
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
    """
    curl: Download a file using curl
    Inputs: url: The URL to retrieve
            cookieURL : The URL of a login page
            username: If logging in, the username
            password: If logging in, the password
            ext: Optionally specify a file extension
    Output: Acquire object
    """
    ext = self.__getExtension(url) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ url, cookieURL, username, password ])
    self.__fileName = self.__dlDir + '/' + curlHash +  ext

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

  def _ftp(self, server, location, username=None, password=None, ext=None):
    """
    ftp: Download a file with FTP
    Inputs: server: The server to access
            location: The path to the file on the server
            username: The username to access the server (if None, default is used)
            password: The password to access the erver
            ext: Optionally specify a file extension
    Output: Acquire object
    """
    # Currently not working...
    ext = self.__getExtension(location) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ server, location, username, password ])
    self.__fileName = self.__dlDir + '/' + curlHash

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    fs.mkdirname(self.__fileName)

    conn = ftplib.FTP(server)
    if (username is not None) and (password is not None):
      conn.login(username, password)
    else:
      conn.login()
    #fi
    p = conn.retrbinary(location, open(self.__fileName, 'wb').write)
    p = int(result.split(' ')[0])
    if p == 226:
      return 0
    #fi
    return p
  #edef

  def _lftp(self, server, location, username, password, ext=None):
    """
    lft: Download a file from an ftp server (but for example with sftp access)
    Inputs: server: The server to access
            location: The path to the file on the server
            username: The username to access the server (if None, default is used)
            password: The password to access the erver
            ext: Optionally specify a file extension
    Output: Acquire object
    """
    ext = self.__getExtension(location) if ext is None else ('.' + ext)
    curlHash = self.__downloadHash([ server, location, username, password ])
    self.__fileName = self.__dlDir + '/' + curlHash

    if self.__checkExistsTag(self.__fileName) and (not self.__redo):
      return 0
    #fi

    if not exe.exists('lftp'):
      msg.error("'lftp' is not installed. Please install in order to continue.")
      return 1
    #fi

    fs.mkdirname(self.__fileName)
    cmd = "echo -en  'open \"%s\"\\nuser \"%s\" \"%s\"\\ncat \"%s\"' | lftp > '%s'" % (server, username, password, location, self.__fileName)
    p = exe.runCommand(cmd, shell=True, verbose=True)
    return p
  #edef


  def _wget(self, url, ext=None):
    """
    wget: Download a file with wget
    Inputs: url: URL to retrieve
            ext: Optionally specify a file extension
    Output: Acquire object
    """
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
    """
    local: Use a locally sourced file
    Inputs: fileName: URI of the local file
    Output: Acquire object
    """
    self.__fileName = fileName
    if os.path.isfile(self.__fileName):
      return 0
    else:
      return 1
    #fi
  #edef

  def _touch(self, fileName):
    """
    touch: Create a local, empty file at a specific location
    Inputs: fileName: URI of the local file
    Output: Acquire object
    """
    self.__fileName = fileName
    return fs.touchFile(fileName)
  #edef

  def _merge(self, acquireObjects, method='cat'):
    """
    merge: Merge multiple acquire objects.
    Inputs: acquireObjects: A list of Acquire objects
            method: How to merge them (implemented: cat, zcat)
    Output: Acquire object
    """
    fileNames = [ ao.acquire() for ao in acquireObjects ]
    if None in fileNames:
      return 1
    #fi

    curlHash = self.__downloadHash(fileNames)
    self.__fileName = self.__dlDir + '/' + curlHash + '.' + method

    if method == 'cat':
      cmd = "cat '%s' > '%s'" % ("' '".join(fileNames), self.__fileName)
    elif method == 'zcat':
      cmd = "cat '%s' | zcat > '%s'" % ("' '".join(fileNames), self.__fileName)
    else:
      raise NotImplementedError("Method '%s' is not implemented for merge" % method)
    #fi
    p = exe.runCommand(cmd, verbose=True, shell=True)
    return p
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

  def _func(self, function):
    oldFile = self.__fileName
    self.__fileName = self.__fileName + '.func'
    return function(oldFile, self.__fileName)
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

    fs.mkdirp(outDirName)

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
      if settings.platform() == "OSX":
        # The -T option doesn't work on Mac
        p = exe.runCommand("cp -R '%s' '%s'" % (oldFile, self.__fileName), verbose=True)
      else:
        p = exe.runCommand("cp -R -T '%s' '%s'" % (oldFile, self.__fileName), verbose=True)
      #fi
    #fi
    return p
  #edef

#eclass
 
