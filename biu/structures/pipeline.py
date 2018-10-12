import hashlib
import json
import tempfile
import json
import os
import threading

from ..config import settings as settings
from .. import utils

snakemake = utils.py.loadExternalModule('snakemake', attr='snakemake')

###############################################################################

def defaultSnakemakeOptions():
  return {
   "use_conda" : True,
   "conda_prefix" : settings.getPipelineCondaPrefixDir()
  }
#edef

def snakemakeOptionsToCmdArgs(opts):
    cmdArgs = []
    for k,v in opts.items():
        if (k == 'use_conda'):
            if v is True:
                cmdArgs.append('--use-conda')
            #fi
        elif k == 'conda_prefix':
            cmdArgs.append("--conda-prefix='%s'" % v)
        # If there are things we don't know, try formatting them with these last two rules.
        elif v is True:
            cmdArgs.append(' --%s' % k.replace('_', '-'))
        else:
            cmdArgs.append(" --%s='%s'" % (k.replace('_', '-'),v))
        #fi
    #efor
    return ' '.join(cmdArgs)
#edef

###############################################################################

def defaultSnakemakeConfig(name):
  return {
    "common_dir" : "%s/%s" % (settings.getPipelineCommonDir(), name),
    "tmp_dir" : "%s/%s" % (settings.getPipelineTemporaryInputDir(), name),
    "biu_settings" : settings.settings,
    "biu_location" : settings.biuLocation,
    "outdir" : None
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
    """
    Run the pipeline.
    Inputs:
    """
    self.__writeConfigFile()
    if threading.current_thread() == threading.main_thread():
        self.__success = snakemake(snakefile=self.__snakefile, configfile=self.__configFile, targets=targets, **self.__snakemakeOptions)
    else:
        utils.msg.warning("There was an error running snakemake through the normal python interface.")
        utils.msg.warning("I will run it via the command line.")
        utils.msg.warning("There may be some problems with translating the python snakemake options to the command line versions.")

        p = utils.exe.runCommand("snakemake --snakefile '%s' --configfile '%s' %s %s"  % (self.__snakefile, self.__configFile, snakemakeOptionsToCmdArgs(self.__snakemakeOptions), ' '.join(targets)))
        self.__success = True if (p == 0) else False
    #etry
  #edef

  def __configHash(self):
    c = sorted([ (k, self.__config[k]) for k in self.__config if k not in defaultSnakemakeConfig("").keys() ], key=lambda x: x[0])
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
