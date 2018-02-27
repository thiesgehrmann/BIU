from ..structures import fileManager as fm
from ..structures import resourceManager as rm

###############################################################################

versions = { "human" : {
  "gff" : "https://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Homo%20sapiens%20(Human)%20\[9606\]%22&format=gff&force=yes"
  },
  "mouse" : {
  "gff" : "https://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Mus%20musculus%20(Mouse)%20\[10090\]%22&format=gff&force=yes"
  }
}

def urlFileIndex(version):
  files = {}
  files["gff"] = (versions[version]["gff"], 'uniprot_%s.gff3' % version, {"urlIsGzipped":True})
  return { k : (u, 'uniprot_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class UniProt(fm.FileManager):

  version = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["annots"], **kwargs)
    self.version = version

    self.annots = rm.GFF3ResourceManager(self, "gff",
                                         parentField=lambda e: e.seqid, 
                                         idField=lambda e: e.attr["ID"] if "ID" in e.attr else '%s;%s;%d;%d' % (e.seqid, e.feature, e.start, e.end),
                                         allowAdditionalColumns=True)

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  #############################################################################

  def getProteinDomains(self, protein):
    return self.annots.getChildren(protein, [ "Domain", "domain"])


#eclass
