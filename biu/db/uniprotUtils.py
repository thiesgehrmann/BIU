from ..structures import Dataset
from .. import formats
from .. import utils

from ..config import settings

###############################################################################

class UniProt(Dataset):

  version = None

  versions = {
    "human" : {
      "gff" : "https://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Homo%20sapiens%20(Human)%20\[9606\]%22&format=gff&force=yes"
    },
    "mouse" : {
      "gff" : "https://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Mus%20musculus%20(Mouse)%20\[10090\]%22&format=gff&force=yes"
    }
  }

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fileIndex = self.__genFileIndex(version)
    Dataset.__init__(self, fileIndex, **kwargs)
    self.version = version

    self._registerObject('_annots', formats.GFF3, ['gff'], fileIndex['gff'].path,
                                         parentField=lambda e: e.seqid,
                                         idField=lambda e: e.attr["ID"] if "ID" in e.attr else '%s;%s;%d;%d' % (e.seqid, e.feature, e.start, e.end),
                                         allowAdditionalColumns=True)

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/uniprot/%s' % ( (settings.getDataDir() if where is None else where), version)
     vData = self.versions[version]
     files = {}
     files['gff'] = utils.Acquire(where=where).curl(vData['gff']).gunzip().finalize('%s/gff.gff3' % finalPath)
     return files
  #edef

  #############################################################################

  def getProteinDomains(self, protein):
    return self._annots.getChildren(protein, [ "Domain", "domain"])
  #edef


#eclass
