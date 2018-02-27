import json
import inspect, os

class globalSettings(object):

  __FILE_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

  _settingsFile = '%s/config.json' % __FILE_DIR__

  _settings = None

###############################################################################

  def __init__(self):
    self._settings = json.load(open(self._settingsFile, "r"))

  def loadSettings(self, fileName):
    self._settings.update(json.load(open(fileName, "r")))
  #edef
  
  def setSettings(self, **kwargs):
    self._settings.update(kwargs)
  #edef
  
  def dumps(self):
    return json.dumps(self._settings)
  #edef
  
  def getSetting(self, settingID):
    if settingID not in self._settings:
      return None
    else:
      return self._settings[settingID]
    #fi
  #edef

  ###############################################################################

  def getWhere(self):
    return self.getSetting("where")
  #edef

  def setWhere(self, where):
    self.setSettings(where=where)
  #edef
  
  ###############################################################################
  
  def getNeo4jDir(self):
    return self.getSetting("neo4j_install_dir")
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
    return self.getSettings("debug_stream")
  #edef

  def setDebugStream(self, stream):
    if stream not in [ "stderr", "stdout" ]:
      stream = "stdout"
    #fi
    return self.setSettings(debug_stream=stream)
  #edef

  def __str__(self):
    dstr = "Configuration:\n"
    for settingID in self._settings:
      dstr += " %s : '%s'\n" % (settingID, self.getSetting(settingID))
    #efor
    return dstr
  #edef

#eclass

settings = globalSettings()
