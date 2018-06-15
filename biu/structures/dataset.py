import os

from .. import utils

class Dataset(object):

  def __init__(self, fileIndex, localCopy={}, overwrite=False):

    self.__fileIndex = fileIndex
    for what in localCopy:
      fileName = localCopy[what]
      self.__fileIndex[what] = utils.Acquire(fileName)
    #efor

    if overwrite:
      for what in self.__fileIndex:
        self.__fileIndex[what] = self.__fileIndex[what].redo()
      #efor
    #fi

    self.__loadedObjects = {}
    self.__registeredObjects = {}
    self.__str_functions = []
  #edef

  def _registerObject(self, oname, objectFormat, requiredFiles, *pargs, **kwargs):
    self.__registeredObjects[oname] = ( objectFormat, requiredFiles, pargs, kwargs )
  #edef

  def _getObject(self, oname):
    if oname in self.__loadedObjects:
      return self.__loadedObjects[oname]
    else:
      return self.__loadObject(oname)
    #fi
  #edef

  def _objectExists(self, oname):
    return oname in self.__registeredObjects
  #edef

  def __loadObject(self, oname):
    if oname not in self.__registeredObjects:
      utils.msg.error("Object '%s' not in dataset." % oname)
      return None
    #fi

    # Acquire the requisite files
    oformat, requiredFiles, pargs, kwargs = self.__registeredObjects[oname]
    for what in requiredFiles:
      if self.__fileIndex[what].acquire() is None:
        utils.msg.error("Could not load object '%s', because requirement '%s' failed." % (oname, what))
        return None
      #fi
    #efor

    # If that was succesful, load the object
    self.__loadedObjects[oname] = oformat(*pargs, **kwargs)
    return self.__loadedObjects[oname]
  #edef

  #############################################################################

  def _addStrFunction(self, func):
    self.__str_functions.append(func)
  #edef

  #############################################################################

  def __str__(self):
    dstr  = "%s object\n" % self.__class__.__name__

    for f in self.__str_functions:
      fstr = f(self)
      if fstr is None:
        continue
      #fi
      for line in fstr.split('\n'):
        dstr += ' ' + line + '\n'
    #efor

    dstr += ' Objects:\n'
    for oname in self.__registeredObjects:
      loaded = oname in self.__loadedObjects
      dstr += '  * [%s] %s\n' % (('X' if loaded else ' '), oname)
    #efor

    dstr += " Files:\n"
    for what in self.__fileIndex:
      loc = self.__fileIndex[what].path
      if os.path.islink(loc):
        dstr += "  * [%s] %s : %s -> %s\n" % ('S' if self.__fileIndex[what].exists else ' ', what, loc, Path(loc).resolve())
      else:
        dstr += "  * [%s] %s : %s\n" % ('X' if self.__fileIndex[what].exists else ' ', what, loc)
      #fi
    #efor
    return dstr
  #edef 

  def __getattr__(self, oname):
    if not self._objectExists(oname):
      raise NameError(oname)
    #fi
    return self._getObject(oname)
  #edef

  #############################################################################

#############################################################################
