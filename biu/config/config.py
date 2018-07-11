import json
import inspect, os

class globalSettings(object):

  __FILE_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

  __settingsFile = '%s/config.json' % __FILE_DIR__

  __settings = None

  __missingDependencies = None

###############################################################################

  def __init__(self):
    self.__settings = json.load(open(self.__settingsFile, "r"))
    self.__missingDependencies = []
  #edef

  def loadSettings(self, fileName):
    self.__settings.update(json.load(open(fileName, "r")))
  #edef
  
  def setSettings(self, **kwargs):
    self.__settings.update(kwargs)
  #edef
  
  def dumps(self):
    return json.dumps(self.__settings)
  #edef
  
  def getSetting(self, settingID):
    if settingID not in self.__settings:
      return None
    else:
      return self.__settings[settingID]
    #fi
  #edef

  @property
  def settings(self):
    return self.__settings
  #edef

  @property
  def biuLocation(self):
    return os.path.dirname(os.path.abspath(self.__FILE_DIR__ + '../../'))
  #edef

  ###############################################################################

  def getWhere(self):
    return os.path.abspath(self.getSetting("where"))
  #edef

  def setWhere(self, where):
    self.setSettings(where=where)
  #edef

  ###############################################################################

  def getDownloadDir(self):
    dldir = self.getSetting("download_where")
    if self.getSetting("download_where") is '':
      return self.getWhere() + '/_downloads'
    else:
      return os.path.abspath(dldir)
    #fi
  #edef

  def setDownloadDir(self, dirName):
    self.setSettings(download_where=dirName)
  #edef

  ###############################################################################

  def getPipelineOutdir(self):
    path = self.getSetting("pipeline_where")
    if path == '':
      path = self.getWhere()
    #fi
    return '%s/%s' % (os.path.abspath(path), self.getSetting('pipeline_base'))
  #edef

  def setPipelineOutdir(self, outdir):
    self.setSettings(pipelines_outdir_base=outdir)
  #edef

  def getPipelineTemporaryInputDir(self):
    return '%s/%s' % (self.getPipelineOutdir(), self.getSetting("pipelines_temporary_indir_name"))
  #edef

  def getPipelineCondaPrefixDir(self):
    return '%s/%s' % (self.getPipelineOutdir(), self.getSetting("pipelines_conda_prefix_name"))
  #edef

  def getPipelineCommonDir(self):
    return '%s/%s' % (self.getPipelineOutdir(), self.getSetting("pipelines_common_name"))
  #edef

  ###############################################################################
  
  def getDebugState(self):
    return self.getSetting("debug_messages")
  #edef

  def setDebugState(self, state):
    if isinstance(state, bool):
      self.setSettings(debug_messages=state)
    #fi
  #edef

  def toggleDebug(self):
    self.setSettings(debug_messages=not(self.getDebugState()))
  #edef

  def getDebugStream(self):
    return self.getSetting("debug_stream")
  #edef

  def setDebugStream(self, stream):
    if stream not in [ "stderr", "stdout" ]:
      stream = "stdout"
    #fi
    return self.setSettings(debug_stream=stream)
  #edef

  ###############################################################################

  def registerMissingDependency(self, module):
    self.__missingDependencies.append(module)
  #edef
    
  def missingDependencies(self):
    topLevelModules = set([ m.split('.')[0] for m in self.__missingDependencies ])
    required = set(self.getSetting("dependencies")) & topLevelModules
    optional = set(self.getSetting("dependencies_optional")) & topLevelModules
    other    = topLevelModules - (required | optional)

    return required, optional, other
  #edef

  ###############################################################################

  def __str__(self):
    dstr = "Configuration:\n"
    for settingID in self.__settings:
      dstr += " %s : '%s'\n" % (settingID, self.getSetting(settingID))
    #efor
    return dstr
  #edef

#eclass

settings = globalSettings()
