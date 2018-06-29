import json
import inspect, os

class globalSettings(object):

  __FILE_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

  __settingsFile = '%s/config.json' % __FILE_DIR__

  __settings = None

###############################################################################

  def __init__(self):
    self.__settings = json.load(open(self.__settingsFile, "r"))

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
    dldir = self.getSetting("download_dir")
    if self.getSetting("download_dir") is None:
      return self.getWhere() + '/_downloads'
    else:
      return os.path.abspath(dldir)
    #fi
  #edef

  def setDownloadDir(self, dirName):
    self.setSetting(download_dir=dirName)
  #edef

  ###############################################################################

  def getPipelineOutdir(self):
    return os.path.abspath(self.getSetting("pipelines_base"))
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
  
  def getNeo4jDir(self):
    return os.path.abspath(self.getSetting("neo4j_install_dir"))
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

  def __str__(self):
    dstr = "Configuration:\n"
    for settingID in self.__settings:
      dstr += " %s : '%s'\n" % (settingID, self.getSetting(settingID))
    #efor
    return dstr
  #edef

#eclass

settings = globalSettings()
