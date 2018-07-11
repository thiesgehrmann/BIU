from .. import utils

sqlite3 = utils.py.loadExternalModule('sqlite3')
pd      = utils.py.loadExternalModule("pandas")

class SQLite(object):
  _diskloc = None
  _connection = None

  def __init__(self, fileName):
    self._diskloc = fileName
    self.connect()
  #edef

  def connect(self):
    try:
      self._connection = sqlite3.connect(self._diskloc)
      self._cursor     = self._connection.cursor()
    except Exception as e:
      utils.error(d)
    #etry
  #edef

  def execute(self, query, values=None):
    r = None
    if values is None:
      r = self._cursor.execute(query)
    else:
      r = self._cursor.execute(query, values)
    #fi
    self._connection.commit()
    return r
  #edef

  def executemany(self, query, values):
    r = self._cursor.executemany(query, values)
    return r
  #edef

  def getTableNames(self):
    res = self.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tNames = [ r[0] for r in res ]
    return tNames
  #edef

  def getTableStructure(self, tableName):
    if tableName not in self.getTableNames():
      utils.error("Table '%s' does not exist in SQLite object" % tableName)
      return ""
    #fi
    res = self.execute("SELECT sql FROM sqlite_master WHERE name=?;", [tableName])
    return [ r[0] for r in res][0]
  #edef

  def dropTable(self, tableName):
    if tableName not in self.getTableNames():
      utils.error("Table '%s' not present in SQLite object" % tableName)
      return False
    else:
      self.execute("DROP TABLE %s;" % tableName )
    #fi
  #edef

  def createTableFromFile(self, fileName, tableName, **kwargs):
    """
      fileName: should be a file
      tableName: Name of table
      tableFields: Array of field definitions (string)
      primaryKey: Primary key for table
    """

    for chunk in pd.read_csv(fileName, chunksize=5000, **kwargs):
      chunk.to_sql(tableName, self._connection, if_exists='append')
    #efor

  #edef

  def createIndex(self, tableName, fieldName, unique=False):
    self.execute("CREATE %s INDEX %s ON %s(%s);" % ("UNIQUE" if unique else "", "idx_%s_%s" % (tableName, fieldName), tableName, fieldName))
  #edef

  def dropIndex(self, tableName, fieldName):
    self.execute("DROP INDEX IF EXISTS idx_%s_%s;" % (tableName, fieldName))
  #edef

  def __str__(self):
    dstr  = "SQLite Object\n"
    dstr += " Where: %s\n" % self._diskloc
    dstr += " Available: %s\n" % ("Yes" if self._cursor is not None else "No")
    if self._cursor is not None:
      dstr += " Tables:\n"
      for res in self.getTableNames():
        dstr += "  * %s\n" % res
      #efor
    #fi
    return dstr
  #edef

#eclass
