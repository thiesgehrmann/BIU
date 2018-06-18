import pandas as pd
import tabix

from .. import utils

from collections import namedtuple

class Tabix(object):

  def __init__(self, fileName, fieldNames=None):
    self.__fieldNames = fieldNames
    self.__fileName   = fileName
    self.__resource   = tabix.open(fileName)
    if self.__fieldNames is not None:
      self.__namedtuple = namedtuple("tabixTsv", fieldNames)
    #fi
  #edef

  def __str__(self):
    dstr =  "Tabix object\n"
    dstr += " Where: %s\n" % self.__fileName
    if self.__fieldNames is not None:
      dstr += " Fields:\n"
      dstr += " * %s" % '\n * '.join(self.__fieldNames)
    #fi
    return dstr
  #edef

  def __safeTabixWrapper(self, chrom, start, end):
    try:
      return self.__resource.query(str(chrom), int(start), int(end))
    except Exception as e:
      utils.msg.error(e)
      return []
    #etry
  #edef 

  def query(self, seqid, start, end, **kwargs):
    return self.queryRegions([(seqid, start, end)], **kwargs)
  #edef

  def queryRegions(self, regions, pandas=False, namedtuple=False):
    res = [ elem for (seq, start, end) in regions for elem in self.__safeTabixWrapper(seq, start, end) ]
    if pandas:
      return pd.DataFrame(res)
    elif namedtuple and (self.__fieldNames is not None):
      return [ self.__namedtuple(*r[:len(self.__fieldNames)]) for r in res ]
    else:
      return res
    #fi
  #edef
#eclass
