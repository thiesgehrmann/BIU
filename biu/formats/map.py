from . import TSVIndex

import pandas as pd

class Map(object):

  __slots__ = [ '__idxs', '__fileName', '__names', '__nameIndex', '__header', '__kwargs', '__tbl' ]

  def __init__(self, fileName, names=None, header=True, **kwargs):
    origNames = names
    if (origNames is None) and header:
      with open(fileName, 'r') as ifd:
        line = ifd.readline()
        names = line.strip().split(kwargs.get('delimiter', '\t'))
      #ewith
    #fi
    self.__header    = header and (origNames is None)
    self.__fileName  = fileName
    self.__names     = names
    self.__nameIndex = dict([ (n, idx) for (idx, n) in enumerate(names)])
    self.__idxs      = {}
    self.__kwargs    = kwargs
    self.__tbl       = None
  #edef

  @property
  def table(self):
    if self.__tbl is None:
      if len(self.__idxs) == 0:
        self.__tbl = self.__getattr__(self.__names[0]).table
      else:
        self.__tbl = self.__idxs.values()[0].table
      #fi
    #fi
    return self.__tbl
  #edef

  def __str__(self):
    dstr = "Indexed TSV Object\n"
    dstr += " Filename: %s\n" % self.__fileName
    dstr += " Indexes:\n"
    for name in self.__names:
      dstr += '  * [%s] %s\n' % ('X' if name in self.__idxs else ' ', name)
    #efor
    
    return dstr
  #edef

  def __getattr__(self, oname):
    if oname not in self.__nameIndex:
      raise NameError(oname)
    #fi
    if oname not in self.__idxs:
      self.__idxs[oname] = TSVIndex(self.__fileName, key=self.__nameIndex[oname], names=self.__names, header=self.__header, **self.__kwargs)  
    #fi
    return self.__idxs[oname]
  #edef

#eclass

