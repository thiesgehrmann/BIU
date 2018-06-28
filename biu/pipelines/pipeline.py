import hashlib
import json
import tempfile
import json
import os

from snakemake import snakemake

from ..config import settings as settings
from .. import utils

###############################################################################

def defaultSnakemakeOptions():
  return {
   "use_conda" : True,
   "conda_prefix" : settings.getPipelineCondaPrefixDir()
  }
#edef

###############################################################################

def defaultSnakemakeConfig(name):
  return {
    "common_dir" : "%s/%s" % (settings.getPipelineCommonDir(), name),
    "tmp_dir" : "%s/%s" % (settings.getPipelineTemporaryInputDir(), name)
  }
#edef

###############################################################################

class Pipeline(object):
  
  __config     = None
  __snakefile  = None
  __success    = False
  __configFile = None
  __autorun    = True
  __snakemakeOptions = None

  def __init__(self, snakefile, config={}, rewriteHashedInputFiles=False, autorun=True, **snakemakeOptions):
    self.__snakefile = snakefile
    self._rewriteHashedInputFiles = rewriteHashedInputFiles
    self.__configFile = None
    self.__autorun = autorun
    self.__snakemakeOptions = defaultSnakemakeOptions().copy()
    self.__snakemakeOptions.update(snakemakeOptions)
    self.__config = defaultSnakemakeConfig(type(self).__name__)
    self.setConfig(config)
  #edef

  def setConfig(self, config=None, **kwargs):
    if config is not None:
      self.__config.update(config)
    #fi
    for (k,v) in kwargs.items():
      self.__config[k] = v
    #fi 

    # Update the hash of the current config
    hashString = self.__configHash()
    self.__config["hash"] = hashString
    self.__config["outdir"] = "%s/%s/%s" % (settings.getPipelineOutdir(), type(self).__name__, hashString)
  #edef

  def __writeConfigFile(self):
    self.__configFile = self.__config["outdir"] + "/config.json"
    utils.mkdirname(self.__configFile)
    with open(self.__config["outdir"] + "/config.json", "w") as ofd:
      json.dump(self.__config, ofd, indent=2)
    #ewith
  #edef

  def run(self, targets=["output"]):
    self.__writeConfigFile()
    self.__success = snakemake(snakefile=self.__snakefile, configfile=self.__configFile, targets=targets, **self.__snakemakeOptions)
  #edef

  def __configHash(self):
    c = sorted([ (k, self.__config[k]) for k in self.__config if k not in ['outdir', 'hash'] ], key=lambda x: x[0])
    c = json.dumps(c)
    return hashlib.md5(c.encode()).hexdigest()
  #edef

  def _generateInputFileName(self, data=None):
    basedir = '%s/%s' % (settings.getPipelineTemporaryInputDir(), type(self).__name__)
    utils.mkdirp(basedir)

    if data is None:
      (_, fileName) = tempfile.mkstemp(dir=basedir)
      exists = False
    else:
      digest = utils.hashArray(data)
      fileName = "%s/digest.%s" % (basedir, digest)
      exists = os.path.isfile(fileName)
    #fi

    if self._rewriteHashedInputFiles:
      exists = False
    #fi

    return fileName, exists
  #edef

  @property
  def success(self):
    return self.__success
  #edef

  @property
  def autorun(self):
    return self.__autorun
  #edef

  @property
  def config(self):
    return self.__config

  def __str__(self):
    dstr  = "%s Pipeline object\n" % type(self).__name__
    dstr += " Where: %s\n" % self.__snakefile
    dstr += " Output: %s\n" % self.__config["outdir"]
    dstr += " Config: %s/config.json\n" % self.__config["outdir"]
    dstr += " Hash: %s\n" % self.__config["hash"]
    dstr += " Completed: %s\n" % ("Yes" if self.__success else "No")
    return dstr
  #edef
#eclass
