import hashlib
import json
import tempfile
import json
import os

from snakemake import snakemake

from ..config import settings as settings
from .. import utils

def defaultSnakemakeOptions():
  return {
   "use_conda" : True,
   "conda_prefix" : settings.getPipelineCondaPrefixDir()
  }
#edef

def defaultSnakemakeConfig(name):
  return {
    "common_dir" : "%s/%s" % (settings.getPipelineCommonDir(), name),
    "tmp_dir" : "%s/%s" % (settings.getPipelineTemporaryInputDir(), name)
  }
#edef

class Pipeline(object):
  
  wf = None # Workflow
  config = None
  snakefile = None
  success = False
  configFile = None

  def __init__(self, snakefile, rewriteHashedInputFiles=False, **kwargs):
    self.snakefile = snakefile
    self._rewriteHashedInputFiles = rewriteHashedInputFiles
  #edef

  def setConfig(self, config=None, configFile=None, **kwargs):

    if configFile is not None:
      self.configFile = configFile
      with open(configFile, "r") as ifd:
        self.config = json.load(ifd)
      #ewith
    elif config is not None:
      self.config = defaultSnakemakeConfig(type(self).__name__)
      self.config.update(config)

      if "outdir" not in self.config:
        hashString = self._configHash()
        self.config["hash"] = hashString
        self.config["outdir"] = "%s/%s/%s" % (settings.getPipelineOutdir(), type(self).__name__, hashString)
      #fi

      self.configFile = self.config["outdir"] + "/config.json"
      utils.mkdirname(self.configFile)
      with open(self.config["outdir"] + "/config.json", "w") as ofd:
        json.dump(self.config, ofd, indent=2)
      #ewith
    else:
      utils.error("You must specify either a config or a configFile.")
      return None
    #fi
  #edef

  def run(self, targets=["output"], **kwargs):
    snakemakeOptions = defaultSnakemakeOptions()
    snakemakeOptions.update(kwargs)

    #smCommand = "snakemake --config '%s' %s" % (self.configFile, ' '.join([ '--%s' % k if 
    self.success = snakemake(snakefile=self.snakefile, config=self.config, targets=targets, **snakemakeOptions)
  #edef

  def _configHash(self):
    c = sorted([ (k, self.config[k]) for k in self.config if k not in ['outdir'] ], key=lambda x: x[0])
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

  def __str__(self):
    dstr  = "Pipeline object for %s\n" % type(self).__name__
    dstr += " Where: %s\n" % self.snakefile
    dstr += " Output: %s\n" % self.config["outdir"]
    dstr += " Config: %s/config.json\n" % self.config["outdir"]
    return dstr
  #edef
#eclass
