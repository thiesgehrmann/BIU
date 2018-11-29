import json
import inspect, os
import sys

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

  def platform(self):
    """
      Return the platform type of the current system [LINUX|OSX|WINDOWS|OTHER]
    """
    platform = sys.platform
    if platform == "linux" or platform == "linux2":
        return "LINUX"
    elif platform == "darwin":
        return "OSX"
    elif platform == "win32":
        return "WINDOWS"
    #fi
    return "OTHER"
  #Edef


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

  def getDataDir(self):
    path = self.getSetting("data_where")
    if path == '':
      path = self.getWhere()
    #fi
    return '%s/%s' % (os.path.abspath(path), self.getSetting('data_base')) 
  #edef

  def setDataDir(self, path):
    self.setSettings(data_where=path)
  #edef

  ###############################################################################

  def getDownloadDir(self):
    path = self.getSetting("download_where")
    if path == '':
      path = self.getWhere()
    #fi
    return '%s/%s' % (os.path.abspath(path), self.getSetting('download_base'))
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
    return '%s/%s' % (os.path.abspath(path), self.getSetting('pipelines_base'))
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
    """Get the state of whether or not debug messages should be displayed"""
    return self.getSetting("debug_messages")
  #edef

  def setDebugState(self, state):
    """Set the state of whether or not debug messages should be displayed"""
    if isinstance(state, bool):
      self.setSettings(debug_messages=state)
    #fi
  #edef

  def toggleDebug(self):
    """Toggle the state of whether or not debug messages should be displayed"""
    self.setSettings(debug_messages=not(self.getDebugState()))
  #edef

  def getDebugStream(self):
    """ Get the stream to which debug messages should be printed to"""
    return self.getSetting("debug_stream")
  #edef

  def setDebugStream(self, stream):
    """ Set the stream to which debug messages should be written to (must be stderr or stdout, falls back to stdout)"""
    if stream not in [ "stderr", "stdout" ]:
      stream = "stdout"
    #fi
    return self.setSettings(debug_stream=stream)
  #edef

  ###############################################################################

  def getErrorState(self):
    """Get the state of whether or not error messages should be displayed"""
    return self.getSetting("error_messages")
  #edef
  def setErrorState(self, state):
    """Set the state of whether or not error messages should be displayed"""
    if isinstance(state, bool):
      self.setSettings(error_messages=state)
  #edef
  def toggleErrorState(self):
    """Toggle the state of whether or not error messages should be displayed"""
    self.setSettings(error_messages=not(self.getErrorState()))

  def getWarningState(self):
    """Get the state of whether or not warning messages should be displayed"""
    return self.getSetting("warning_messages")
  #edef
  def setWarningState(self, state):
    """Set the state of whether or not warning messages should be displayed"""
    if isinstance(state, bool):
      self.setSettings(warning_messages=state)
  #edef
  def toggleWarningState(self):
    """Toggle the state of whether or not warning messages should be displayed"""
    self.setSettings(warning_messages=not(self.getWarningState()))

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
