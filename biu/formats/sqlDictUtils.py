from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

import os
import json

def urlFileIndex(name):
  files = {}

  files["sqlite_db"] = (None, "sqlDict/%s.sqlDict.sqlite" % name, {})

  return files
#edef

class SQLDict(fm.FileManager):
  """SQLDict is designed to behave like a dictionary, except that it stores the values of the dictionary in JSON Strings in a SQLite database behind the scenes.
  To improve speed, it also caches the values that it stores and retrieves from the SQLite database during runtime.
  This allows you to keep a running dataset of values without the need to recompute stuff each time.

  Operations supported:

  Initialization:
    D = SQLite("mydict")
    D["key"] = {"a" : [ 1, 2, 3, "faf", {1: 2}], "5" : -1 }
    for key in D:
      print(D[key])
    if "key" in D:
      print(D["key"])
  """

  _sqlDict = None
  _cache   = None

  def __init__(self, name, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(name), objects=[ "_sqlDict" ], **kwargs)
    self._sqlDict = rm.SQLiteResourceManager(self, "sqlite_db")
    self._cache   = {}

    if not(self.haveFile("sqlite_db")):
      self.touchFile("sqlite_db")
      self._sqlDict.execute("CREATE TABLE data(id STRING PRIMARY KEY, value TEXT);")
    #fi
  #edef

  def _store(self, key, value):
    self._cache[key] = value
    return self._sqlDict.execute("REPLACE INTO data(id, value) VALUES (?, ?);", [str(key), json.dumps(value)])
  #edef

  def _retrieve(self, key):
    key = str(key)

    if key in self._cache:
      return self._cache[key]
    #fi

    res = list(self._sqlDict.execute("SELECT value FROM data WHERE id IS ?;", [key]))
    if len(res) == 0:
      return None
    else:
      res = json.loads(res[0][0])
      self._cache[key] = res
      return res
    #fi
  #edef

  def __getitem__(self, key):
    return self._retrieve(key)
  #edef

  def __setitem__(self, key, value):
    return self._store(key, value)
  #edef

  def __contains__(self, key):
    key = str(key)
    if (key in self._cache) or (self.__getitem__(key) is not None):
      return True
    else:
      return False
    #fi
  #edef
 
  def __iter__(self):
    self._iterKeys = list(self._sqlDict.execute("SELECT id FROM data;"))
    return self
  #edef

  def __next__(self):
    if len(self._iterKeys) == 0:
      raise StopIteration
    #fi
    v = self._iterKeys.pop()
    return v[0]
  #edef
 
#eclass
