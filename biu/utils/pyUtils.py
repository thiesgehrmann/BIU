import types
import importlib
import importlib.machinery
import inspect

from . import msgUtils as msg
from ..config import settings

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

def loadExternalModule(module, attr=None):
  class AbsentModule(object):
    def __init__(self, module, exception):
      self.__module = module
      self.__exception = exception
    #edef

    def __re__(self, *pargs, **kwargs):
      msg.error("For this functionality, you need to install '%s'" % self.__module)
      raise ModuleNotFoundError from self.__exception
    #edef
    
    __getattr__ = __re__
    __getitem__ = __re__
    __call__    = __re__
  #eclass

  lmod = None
  try:
    lmod = importlib.import_module(module)
    if attr is not None:
      lmod = getattr(lmod, attr)
    #fi
  except ModuleNotFoundError as e:
    settings.registerMissingDependency(module)
    lmod = AbsentModule(module, e)
  #etry
  return lmod
#edef

def source(obj):
  if isinstance(obj, str):
    msg.error("Not implemented")
    return None
  else:
    return inspect.getsource(obj)
  #fi
#edef
