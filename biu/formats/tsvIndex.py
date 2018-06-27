from .. import utils

import csv
from collections import namedtuple
import pandas as pd

class TSVIndex(object):

  __slots__ = [ '__idx', '__lst', '__key', '__fileName', '__tbl', '__names' ]

  def __init__(self, fileName, key=0, names=None, **kwargs):
    self.__idx = {}
    self.__lst = []
    self.__fileName = fileName
    self.__tbl = None
    self.__names = names

    if isinstance(key, int):
      self.__key = tuple([key])
    #fi

    if names is not None:
      namedTupleObject = namedtuple('TSVIndexRow', names)
    #fi

    itemKey = lambda item: ','.join([ item[i] for i in self.__key])

    with open(fileName, 'r') as csvfile:
      reader = csv.reader(csvfile, **kwargs)
      for i, row in enumerate(reader):  
         self.__lst.append( namedTupleObject(*row) if names is not None else row )
         self.__idx[itemKey(row)] = self.__idx.get(itemKey(row), []) + [ i ]
      #efor
    #ewith
  #edef

  @property
  def table(self):
    if self.__tbl is None:
      self.__tbl = pd.DataFrame(self.__lst, columns=self.__names)
    #fi
    return self.__tbl
  #edef

  def lookup(self, key, singleton=False):
    if key not in self.__idx:
      utils.msg.error("Item '%s' not in map." % key)
      return None
    #fi

    if singleton:
      return self.__lst[self.__idx[key][0]]
    #fi

    return [ self.__lst[i] for i in self.__idx[key] ]
  #edef

  def __getitem__(self, key):
    return self.lookup(key)
  #edef

  def __iter__(self):
    return self.__lst.__iter__()
  #edef

  def keys(self):
    return self.__idx.keys()
  #edef

  def values(self):
    return self.__lst()
  #edef

  def __str__(self):
    dstr = "Indexed TSV Object\n"
    dstr += " Filename: %s\n" % self.__fileName
    dstr += " Indexed on column %s\n" % (str(self.__key))
    dstr += " #Rows : %d\n" % len(self.__lst)
    dstr += " #Indexes: %d\n" % len(self.__idx)
    return dstr
  #edef

#eclass
          
