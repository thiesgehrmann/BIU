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
  },

  "Ensembl_GRCh37" : {
    "genomeURLs" : { "all" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz" },
    "gffURL" : "ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz",
    "chr" : [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "X", "Y" ],
    "cdsURL" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz" },

  "RefSeq_GRCh37" : {
    "genomeURLs" : { "all" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz" },
    "gffURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz",
    "chr" : [ "all" ],
    "cdsURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_rna.fna.gz",
    "aaURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_protein.faa.gz" },

  "RefSeq_GRCh38" : {
    "genomeURLs" : { "all" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" },
    "gffURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz",
    "chr" : [ "all" ],
    "cdsURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz",
    "aaURL" : "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz" }
          
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

    def loadedObjects(s):
      dstr = 'Objects:\n'
      dstr += " * [%s] _gffObject\n" % ('X' if s._gffObject is not None else ' ')
      return dstr
    #edef

    self.str_functions.append( lambda s: loadedObjects(s))
    self.str_functions.append( lambda s: "Genome : %s" % s.genomeID )
  #edef

  def __urlFileIndex(self):
    files = {}
    if "gffURL" in genomes[self.genomeID]:
      files["gff"] = ( genomes[self.genomeID]["gffURL"], self.where + '/genome.gff3', {})
    #fi
    if "cdsURL" in genomes[self.genomeID]:
      files["cds"] = ( genomes[self.genomeID]["cdsURL"], self.where + '/cds.fa', {})
    #fi
    if "aaURL" in genomes[self.genomeID]:
      files["aa"] = ( genomes[self.genomeID]["aaURL"], self.where + '/aa.fa', {})
    if "genomeURLs" in genomes[self.genomeID]:
      for chrID in genomes[self.genomeID]["genomeURLs"]:
        files["chr_%s" % chrID] = (genomes[self.genomeID]["genomeURLs"][chrID], self.where + '/chr%s.fa.gz' % chrID, {})
      #efor
    #fi
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

  ###############################################################################

  

#eclass


