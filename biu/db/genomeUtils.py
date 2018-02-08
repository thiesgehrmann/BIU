import os
import imp

from . import fileManager as fm
from .. import utils as utils
from .. import gff3Utils as gffUtils
imp.reload(utils)
imp.reload(fm)

genomes = {
  "GRCh37" : {
    "genomeURLs" : { chr : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.%s.fa.gz" % chr for chr in [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "X", "Y" ] }, 
    "gffURL" : "ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz",
    "chr" : [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "X", "Y" ],
    "cdsURL" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz"
  }
}

def listGenomes():
  print("Available genomes:")
  for genome in genomes:
    print(" * %s" % genome)
  #efor
#edef

###############################################################################

class Genome(fm.FileManager):

  genomeID = None
  where    = None
  fileIndex = None

  def __init__(self, genomeID, where='./', **kwargs):
    fm.FileManager.__init__(self, where, **kwargs)
    self.genomeID = genomeID
    self.fileIndex = self.__urlFileIndex()

    self.str_functions.append( lambda s: "Genome : %s" % s.genomeID )
  #edef

  def __urlFileIndex(self):
    files = {}
    files["gff"] = ( genomes[self.genomeID]["gffURL"], self.where + '/genome.gff3', {})
    files["cds"] = ( genomes[self.genomeID]["cdsURL"], self.where + '/genome.fa', {})
    for chrID in genomes[self.genomeID]["chr"]:
      files["chr_%s" % chrID] = (genomes[self.genomeID]["genomeURLs"][chrID], self.where + '/chr%s.fa.gz' % chrID, {})
    #efor
    return files
  #edef

  ###############################################################################

  _gffObject = None

  def _requireGffObject(self):
    if self._gffObject is None:
      if not(self.satisfyRequiredFiles(["gff"])):
        return False
      #fi
      self._gffObject = gffUtils.GFF3(fileName=self.getFileName("gff"))
    #fi
    return True
  #edef

  def getGFF(self):
    if not(self._requireGffObject()):
      return None
    else:
      return self._gffObject
    #fi
  #edef

#eclass


