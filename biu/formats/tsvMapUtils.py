from .. import utils
import os

import pandas as pd
import json

class TSVMap(object):
  __mapping = None
  __mappingR = None
  def __init__(self, tsvFile, mapFrom=0, mapTo=1, pickle=True, overwritePickle=False, **kwargs):
    self.__mapping = {}
    self.__mappingR = {}

    print(tsvFile)
    self.__resource = pd.read_csv(tsvFile, **kwargs)

    mapPickleFile  = tsvFile + '.tsvMap.pkl'
    mapRPickleFile = tsvFile + '.tsvMap.r.pkl'

    if pickle and not(overwritePickle) and os.path.exists(mapPickleFile) and os.path.exists(mapRPickleFile):
      utils.msg.dbm("Loading the index from pickle")
      with open(mapPickleFile, 'r') as ifd:
        self.__mapping  = json.load(ifd)
      with open(mapRPickleFile, 'r') as ifd:
        self.__mappingR = json.load(ifd)
    else:
      utils.msg.dbm("Generating the index")
      for row in self.__resource.itertuples():
        i = row.Index
        fromValue = str(row[mapFrom+1])
        toValue   = str(row[mapTo+1])
        if fromValue not in self.__mapping:
          self.__mapping[fromValue] = []
        #fi
        if toValue not in self.__mappingR:
          self.__mappingR[toValue] = []
        #fi
        self.__mapping[fromValue].append((toValue, i))
        self.__mappingR[toValue].append((fromValue, i))
      #efor
      if pickle:
        utils.msg.dbm("Pickling the index")
        with open(mapPickleFile, 'w') as ofd:
          json.dump(self.__mapping, ofd)
        #ewith
        with open(mapRPickleFile, 'w') as ofd:
          json.dump(self.__mappingR, ofd)
        #ewith
      #fi
    #fi
  #edef

  def __lookup(self, key, inverse, withEntry):
    mapping = self.__mappingR if inverse else self.__mapping
    if key not in mapping:
      utils.msg.error("'%s' not in map" % key)
      return []
    else:
      if withEntry:
        return mapping[key]
      else:
        return [ v[0] for v in mapping[key] ]
      #fi
    #fi
  #edef

  def __getitem__(self, key, withEntry=False):
    return self.__lookup(key, False, withEntry=withEntry)
  #edef

  def lookup(self, key, withEntry=False):
    return self.__lookup(key, False, withEntry=withEntry)

  def inverse(self, key, withEntry=False):
    return self.__lookup(key, True, withEntry=withEntry)
  #edef

  @property
  def fromKeys(self):
    return list(self.__mapping.keys())
  #edef

  @property
  def toKeys(self):
    return list(self.__mappingR.keys())
  #edef

  def invert(self):
    mapping = self._mapping
    self.mapping = self._mappingR
    self.mappingR = mapping
  #edef

#eclass
