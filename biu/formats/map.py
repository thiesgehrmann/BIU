from .. import utils

from . import TSVIndex

pd = utils.py.loadExternalModule("pandas")

class Map(object):
  """
  Index columns in a tsv file.
  Each column is indexed using a TSVIndex formatted object.
  """

  __slots__ = [ '__idxs', '__fileName', '__names', '__nameIndex', '__onlyIndex', '__header', '__kwargs', '__tbl' ]

  def __init__(self, fileName, names=None, onlyIndex=None, header=True, **kwargs):
    """
    Initialize the Map object.
    
    parameters:
    -----------
    fileName: The path to the relevant file
    names:    The names of the fields
    only_index: Only index the columns specified in this argument.
    header:   Is there a header in the file?
    **kwargs: Additional arguments to TSVIndex (e.g. sep, etc).
    
    """
    origNames = names
    if (origNames is None):
      with open(fileName, 'r') as ifd:
        line = ifd.readline()
        lineValues = line.strip().split(kwargs.get('delimiter', '\t'))
      #ewith
      if header:
        names = lineValues
      else:
        names = [ 'f%d' % f for f in range(len(lineValues)) ]
      #fi
    #fi

    if onlyIndex is None:
      self.__onlyIndex = names
    else:
      self.__onlyIndex = list(set(onlyIndex) & set(names))
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
    """
    Return a table of the mappings
    """
    if self.__tbl is None:
      if len(self.__idxs) == 0:
        self.__tbl = getattr(self, self.__onlyIndex[0]).table
      else:
        self.__tbl = list(self.__idxs.values())[0].table
      #fi
    #fi
    return self.__tbl
  #edef

  def __str__(self):
    """
    String representation of object.
    """
    dstr = "Indexed TSV Map Object\n"
    dstr += " Filename: %s\n" % self.__fileName
    dstr += " Indexes:\n"
    for name in self.__onlyIndex:
      dstr += '  * [%s] %s\n' % ('X' if name in self.__idxs else ' ', name)
    #efor
    
    return dstr
  #edef
    
  def __repr__(self):
    """
    String representation of object.
    """
    return str(self)
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
    
  def __dir__(self):
      """
      To allow tab-completion for the registered objects.
      """
      return object.__dir__(self) + list(self.__nameIndex.keys())
  #edef

#eclass

