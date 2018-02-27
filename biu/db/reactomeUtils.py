from ..structures import fileManager as fm
from ..structures import resourceManager as rm

from ..config import settings as settings

###############################################################################

versions = { "current" : {
  "neo4j" : "https://neo4j.com/artifact.php?name=neo4j-community-3.3.3-unix.tar.gz",
  "db" : "https://reactome.org/download/current/reactome.graphdb.tgz",
  "user" : "neo4j",
  "passwd" : "neo4j",
  "default_password" : "neo4j",
  "dbName" : "graph.db"
  }
}

def urlFileIndex(version):
  files = {}
  files["neo4j"] = (versions[version]["neo4j"], 'neo4j.tgz', {})
  files["db"]    = (versions[version]["db"], 'reactome.graphdb.tgz', {}) 
  return { k : (u, 'reactome_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class Reactome(fm.FileManager):

  version = None
  instance = None
  connection = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), ["instance", "connection"], **kwargs)
    self.version = version
    self.instance = rm.Neo4JInstanceResourceManager(self, "neo4j",
                                                          location=self._where)
    self.connection = rm.Neo4JDBResourceManager(self, "db", self.instance, 
                                                username=versions[self.version]["user"], 
                                                password=versions[self.version]["passwd"], 
                                                dbName=versions[self.version]["dbName"])
                                                #changeFromDefaultPassword=versions[self.version]['default_password'])

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def query(self, *pargs, **kwargs):
    return self.connection.query(*pargs, **kwargs)
  #edef

  def getPathway(self, pathwayID):
    return self.query('MATCH (pathway:Pathway{stId:"%s"}) RETURN pathway' % pathwayID)
  #edef

  def getProteinsInPathway(self, pathwayID):
    return self.query("""MATCH (p:Pathway{stId:"%s"})-[:hasEvent*]->(rle:ReactionLikeEvent),
      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
      (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
      RETURN DISTINCT re.identifier AS Identifier, rd.displayName AS Database""" % pathwayID)
  #edef

  #############################################################################


#eclass
