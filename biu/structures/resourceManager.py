
# Python modules
import itertools
from collections import namedtuple
import os
import json

# Internal modules
from .. import formats
from .. import utils
from ..config import settings as settings

from .lazyObject import LazyObject

# External modules
import pandas as pd
import numpy as np
import tabix
import vcf
from neo4j.v1 import GraphDatabase

###############################################################################

class ResourceManager(LazyObject):

  _resource = None
  _initialized = None
  _fmObject = None
  _requiredFiles = None

  def __init__(self, fmObject, requiredFiles, **kwargs):
    utils.dbm("Initializing the %s object NOW" % type(self).__name__)
    self._fmObject = fmObject
    self._requiredFiles = requiredFiles
    if not(fmObject.satisfyRequiredFiles(requiredFiles)):
      utils.error("Could not initialize")
      self._initialized = False
    else:
      self._initialized = True
    #fi
  #edef

  def __call__(self):
    return self._resource
  #edef

  def __str__(self):
    dstr =  "ResourceManager Object:\n"
    dstr += " Initialized: %s\n" % ("Yes" if self._initialized else "No")
    dstr += " Resources managed: \n"
    for fileID in self._requiredFiles:
      dstr += "  * %s -> %s\n" % (fileID, self._fmObject.getFileName(fileID))
    #efor
    return dstr
  #edef

#eclass

###############################################################################

class TSVResourceManager(ResourceManager):
  def __init__(self, fmObject, tsvFile, fieldNames=None, delimiter='\t', skiprows=0, **kwargs):
    ResourceManager.__init__(self, fmObject, [ tsvFile ], **kwargs)
    if self._initialized:
      self._resource = pd.read_csv(self._fmObject.getFileName(tsvFile), delimiter=delimiter, names=fieldNames, skiprows=skiprows, **utils.stripkwargs(kwargs))
    #fi
  #edef

  def __getattr__(self, attr):
    return self._resource.__getattr__(attr)
  #edef

  def __getitem__(self, key):
    return self._resource.__getitem__(key)
  #edef

  def __setitem__(self, key, value):
    return self._resource.__setitem__(key, value)
  #edef
#eclass

###############################################################################

class TabixTSVResourceManager(ResourceManager):

  def __init__(self, fmObject, tsvFile, tabixFile, fieldNames=None, **kwargs):
    ResourceManager.__init__(self, fmObject, [ tsvFile, tabixFile ], **kwargs)
    if self._initialized:
      self._resource = tabix.open(self._fmObject.getFileName(tsvFile))
      if fieldNames is not None:
        self._namedtuple = namedtuple("tabixTsv", fieldNames)
      #fi
      self._fieldNames = fieldNames
    #fi
  #edef

  def query(self, seqid, start, end, pandas=False, namedtuple=False, tabixIter=True):
    res = utils.tabixQueryWrapper(self._resource, seqid, start, end)
    if pandas:
      return pd.DataFrame([ r for r in res])
    elif namedtuple and (self._fieldNames is not None):
      return [ self._namedtuple(*r) for r in res ]
    elif tabixIter:
      return res
    else:
      return [ r for r in res ]
    #fi
  #edef
#eclass

###############################################################################

class VCFResourceManager(ResourceManager):
  def __init__(self, fmObject, vcfFile, tabixFile, fieldNames=None, **kwargs):
    ResourceManager.__init__(self, fmObject, [ vcfFile, tabixFile ], **kwargs)
    if self._initialized:
      self._resource = vcf.Reader(filename=self._fmObject.getFileName(vcfFile), compressed=True, **kwargs)
    #fi
  #edef

  def query(self, seqid, start, end,
    filters=None,
    gtFilters=None,
    sampleFilters=None,
    types=None,
    subTypes=None,
    extract=None, **kwargs):
    """ This function is hacked together from the pyVCF _Record and _Call classes. I hope that their definitions don't change anytime soon.
    """

    res = formats.VCF.query(self._resource, seqid, start, end)

    # Filter out variants and samples
    # To filter out samples, we need to modify the _Call.samples and _sample_indexes variables.
    # In my version of pyvcf, the is_filtered function is not defined for some reason, so I add it here.

    res = formats.VCF.filterType(res, types=types)
    res = formats.VCF.filterSubTypes(res, subTypes=subTypes)
    res = formats.VCF.filter(res, gtFilters=gtFilters)
    res = formats.VCF.filterSamples(res, sampleFilters=sampleFilters)

    res = formats.VCF.extract(res, extract=extract)

    return res
  #edef

  def queryRegions(self, regions, **kwargs):
    R = []
    for (seqid, start, end) in regions:
      for r in self.query(seqid, start, end, **kwargs):
        R.append(r)
    #efor
    return R
  #edef

###############################################################################

class GFF3ResourceManager(ResourceManager, formats.GFF3):
  def __init__(self, fmObject, gff3File, **kwargs):
    ResourceManager.__init__(self, fmObject, [ gff3File ])
    if self._initialized:
      formats.GFF3.__init__(self, fileName=self._fmObject.getFileName(gff3File), **kwargs)
    #fi
  #edef

  def __str__(self):
    return formats.GFF3.__str__(self)
  #edef
#eclass

###############################################################################

class XLSXResourceManager(ResourceManager, formats.XLSX):
  def __init__(self, fmObject, xlsxFile, **kwargs):
    ResourceManager.__init__(self, fmObject, [ xlsxFile ])
    if self._initialized:
      formats.XLSX.__init__(self, fileName=self._fmObject.getFileName(xlsxFile), **kwargs)
    #fi
  #edef

  def __str__(self):
    return formats.XLSX.__str__(self)
  #edef
#eclass

###############################################################################

class FastaResourceManager(ResourceManager, formats.Fasta):
  def __init__(self, fmObject, fastaFile, fileSkipLines=0, fileMaxLines=None, **kwargs):
    ResourceManager.__init__(self, fmObject, [ fastaFile ], **kwargs)
    if self._initialized:
      formats.Fasta.__init__(self, fileName=self._fmObject.getFileName(fastaFile), **kwargs)
    #fi
  #edef

  def __str__(self):
    return formats.Fasta.__str__(self)
  #edef
#eclass

###############################################################################

class SQLiteResourceManager(ResourceManager, formats.SQLite):

  def __init__(self, fmObject, sqliteFile, **kwargs):
    ResourceManager.__init__(self, fmObject, [ sqliteFile ], **kwargs)
    if self._initialized:
      formats.SQLite.__init__(self, self._fmObject.getFileName(sqliteFile))
    #fi
  #edef
    
  def __str__(self):
    return formats.SQLite.__str__(self)
  #edef
#eclass

###############################################################################

class GAFResourceManager(ResourceManager, formats.GAF):
  def __init__(self, fmObject, gafFile, **kwargs):
    ResourceManager.__init__(self, fmObject, [ gafFile ], **kwargs)
    if self._initialized:
      formats.GAF.__init__(self, self._fmObject.getFileName(gafFile), **kwargs)
    #fi
  #edef

  def __str__(self):
    return formats.GAF.__str__(self)
  #edef
#eclass

###############################################################################

class TSVMapResourceManager(TSVResourceManager):
  _mapping = None
  _mappingR = None
  def __init__(self, fmObject, tsvFile, mapFrom=0, mapTo=1, pickle=True, overwritePickle=False, **kwargs):
    TSVResourceManager.__init__(self, fmObject, tsvFile, index_col=False, **kwargs)
    if self._initialized:
      self._mapping = {}
      self._mappingR = {}

      mapPickleFile  = self._fmObject.getFileName(tsvFile) + '.pkl'
      mapRPickleFile = self._fmObject.getFileName(tsvFile) + '.r.pkl'

      if pickle and not(overwritePickle) and os.path.exists(mapPickleFile) and os.path.exists(mapRPickleFile):
        utils.dbm("Loading the index from pickle")
        with open(mapPickleFile, 'r') as ifd:
          self._mapping  = json.load(ifd)
        with open(mapRPickleFile, 'r') as ifd:
          self._mappingR = json.load(ifd)
      else:
        utils.dbm("Generating the index")
        for row in self._resource.itertuples():
          i = row.Index
          fromValue = str(row[mapFrom+1])
          toValue   = str(row[mapTo+1])
          if fromValue not in self._mapping:
            self._mapping[fromValue] = []
          #fi
          if toValue not in self._mappingR:
            self._mappingR[toValue] = []
          #fi
          self._mapping[fromValue].append((toValue, i))
          self._mappingR[toValue].append((fromValue, i))
        #efor
        if pickle:
          utils.dbm("Pickling the index")
          with open(mapPickleFile, 'w') as ofd:
            json.dump(self._mapping, ofd)
          #ewith
          with open(mapRPickleFile, 'w') as ofd:
            json.dump(self._mappingR, ofd)
          #ewith
        #fi
      #fi
    #fi
  #edef

  def _lookup(self, key, inverse, withEntry):
    mapping = self._mappingR if inverse else self._mapping
    if key not in mapping:
      utils.error("'%s' not in map" % key)
      return []
    else:
      if withEntry:
        return mapping[key]
      else:
        return [ v[0] for v in mapping[key] ]
      #fi
    #fi
  #edef

  #def __getitem__(self, key, withEntry=False):
  #  return self._lookup(key, False, withEntry=withEntry)
  #edef

  def lookup(self, key, withEntry=False):
    return self._lookup(key, False, withEntry=withEntry)

  def inverse(self, key, withEntry=False):
    return self._lookup(key, True, withEntry=withEntry)
  #edef

  def lookupKeys(self):
    return list(self._mapping.keys())
  #edef

  def inverseKeys(self):
    return list(self._mappingR.keys())
  #edef

  def invert(self):
    mapping = self._mapping
    self.mapping = self._mappingR
    self.mappingR = mapping
  #edef
      

###############################################################################

class Neo4JInstanceResourceManager(ResourceManager):

  _bindir = None
  _location = None

  _unpackedFile = 'neo4j_unpacked'
  _restoreFile  = 'neo4j_restoreFile'

  def __init__(self, fmObject, neo4jTarGzFile, location=settings.getNeo4jDir(), overwrite=False):
    ResourceManager.__init__(self, fmObject, [ neo4jTarGzFile ])
    self._location = location
    if self._initialized:
      if overwrite:
        utils.runCommand("rm -rf '%s'" % self._location)
      #fi

      if not(os.path.exists(self._location)):
        utils.mkdirp(self._location)
      #fi

      if not(os.path.isdir(self._location)):
        utils.error("'%s' is not a directory" % self._location)
        return None
      #fi

      if not(os.path.exists('%s/%s' % (self._location, self._unpackedFile))):
        utils.dbm("Untarring...")
        utils.runCommand("tar --overwrite -xvf '%s' -C '%s' --strip-components 1" % (self._fmObject.getFileName(neo4jTarGzFile), self._location))
        utils.dbm("writing status files")
        utils.touchFile('%s/%s' % (self._location, self._unpackedFile))
        utils.touchFile('%s/%s' % (self._location, self._restoreFile))
      #fi

      self._bindir = '%s/bin' % self._location
    #fi
  #edef

  def start(self):
    utils.runCommand('"%s/neo4j" start' % self._bindir)
  #edef

  def stop(self):
    utils.runCommand('"%s/neo4j" stop' % self._bindir)
  #edef

  def restore(self, fileName, dbName, username=None, password=None, hostname=None, port=None, protocol=None, changeFromDefaultPassword=None):
    utils.dbm("Loading database from '%s'" % fileName)
    self.stop()
    utils.runCommand("'%s/neo4j-admin' load --from='%s' --database='%s' --force" % (self._bindir, fileName, dbName))
    with open('%s/%s' % (self._location, self._restoreFile), "w") as ofd:
      ofd.write("%s" % fileName)
    #ewith
    if changeFromDefaultPassword is not None:
      utils.dbm("Changing password from '%s' to '%s'" % (changeFromDefaultPassword, password))
      if (username is None) or (password is None) or (hostname is None) or (port is None) or (protocol is None):
        utils.error("Can't change password without connection information")
      else:
        self.start()
        utils.dbm("%s://%s:%s" % (protocol, hostname, port), (username, changeFromDefaultPassword))
        with GraphDatabase.driver("%s://%s:%s" % (protocol, hostname, port), auth=(username, changeFromDefaultPassword)) as driver:
          with driver.session() as session:
            session.run("CALL dbms.changePassword('%s')" % password)
          #ewith
        #ewith
      #fi
  #edef

  def getRestoreFile(self):
    with open('%s/neo4j_restoreFile' % self._location, "r") as ifd:
      return ifd.read().replace('\n', '')
    #ewith
  #edef

  def status(self):
    (out, err) = utils.getCommandOutput('"%s/neo4j" status' % self._bindir)
    return out
  #edef

  def __exit__(self, exc_type, exc_value, traceback):
    self.stop()
  #edef

#eclass
    
###############################################################################

class Neo4JDBResourceManager(ResourceManager):
  _protocol = None
  _hostname = None
  _port = None
  _username = None
  _password = None
  _driver = None
  _dbName = None
  _neo4jResource = None
  _session = None

  def __init__(self, fmObject, dbFile, neo4jRM, username, password, dbName="graph.db", hostname="localhost", port=7687, protocol="bolt", changeFromDefaultPassword=None):
    ResourceManager.__init__(self, fmObject, [ dbFile ])
    if self._initialized:
      self._protocol = str(protocol)
      self._hostname = str(hostname)
      self._port = str(port)
      self._username = username
      self._password = password
      self._dbName = dbName
      self._neo4jResource = neo4jRM

      if self._neo4jResource.getRestoreFile() != self._fmObject.getFileName(dbFile):
        self._neo4jResource.restore(self._fmObject.getFileName(dbFile), self._dbName, self._username, self._password, self._hostname, self._port, self._protocol, changeFromDefaultPassword)
      #fi
      self._neo4jResource.start()
      utils.dbm("%s://%s:%s" % (self._protocol, self._hostname, self._port))
      self.openConnection()
    #fi
  #edef

  def openConnection(self):
    self._driver = GraphDatabase.driver("%s://%s:%s" % (self._protocol, self._hostname, self._port), auth=(self._username, self._password))
    self._session = self._driver.session()
  #edef

  def closeConnection(self):
    self._session.close()
    self._driver.close()
  #edef

  def changePassword(self, newPassword):
    self.query("CALL dbms.changePassword('%s')" % newPassword)
    self._password = newPassword
    self.closeConnection()
    self.openConnection()

  def query(self, queryString):
    return self._session.run(queryString)
  #edef

  def stop(self):
    self._driver.close()
    self._neo4jResource.stop()
  #edef

#  def query(self, q):
    
  #edef
    
 ###############################################################################
