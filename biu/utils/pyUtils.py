import types
import importlib.machinery

###############################################################################

def loadModuleFromFile(fileName, moduleName=None):
  if moduleName is None:
    moduleName = '.'.join(fileName.split('/')[-1].split('.')[:-1])
  #fi
  loader = importlib.machinery.SourceFileLoader(moduleName, fileName)
  mod = types.ModuleType(loader.name)
  loader.exec_module(mod)
  return mod
#edef

###############################################################################
