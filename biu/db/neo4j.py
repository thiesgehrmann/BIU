from ..structures import fileManager as fm
from ..structures import resourceManager as rm

from ..config import settings as settings

###############################################################################

versions = { "3.3.3" : {
  "neo4j" : "https://neo4j.com/artifact.php?name=neo4j-community-3.3.3-unix.tar.gz",
  }
}

def urlFileIndex(version):
  files = {}
  files["neo4j"] = (versions[version]["neo4j"], 'neo4j.tgz', {})
  return files
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class Neo4jDB(fm.FileManager):

  version = None
  instance = None

  def __init__(self, version=list(versions.keys())[0], where=None, **kwargs):
    if where is None:
      where = settings.getNeo4jDir()
    #fi
    fm.FileManager.__init__(self, urlFileIndex(version), where, ["instance"], **kwargs)
    self.version = version
    self.instance = rm.Neo4JInstanceResourceManager(self, "neo4j")

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  #############################################################################


#eclass
