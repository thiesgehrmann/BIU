
from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

versions = {
  "GRCh37" : {
    "genomeURLs" : { chr : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.%s.fa.gz" % chr for chr in [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "MT", "X", "Y" ] }, 
    "gffURL" : "ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz",
    "chr" : [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "MT", "X", "Y" ],
    "cdsURL" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz"
  },

  "Ensembl_GRCh37" : {
    "genomeURLs" : { "all" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz" },
    "gffURL" : "ftp://ftp.ensembl.org/pub/grch37/update/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz",
    "chr" : [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "X", "Y" ],
    "cdsURL" : "ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz" },

  "Ensembl_GRCh38_91" : {
    "cdsURL" : "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",
    "aaURL" :  "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
    "genomeURLs" : { chr : "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.%s.fa.gz" % chr for chr in [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "MT", "X", "Y" ] },
    "chr" : [ "1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22", "MT", "X", "Y" ],
    "gffURL" : "ftp://ftp.ensembl.org/pub/release-91/gff3/homo_sapiens/Homo_sapiens.GRCh38.91.gff3.gz" },

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

def urlFileIndex(version):
  files = {}
  if "gffURL" in versions[version]:
    files["gff"] = ( versions[version]["gffURL"], 'genome.gff3', {})
  #fi
  if "cdsURL" in versions[version]:
    files["cds"] = ( versions[version]["cdsURL"], 'cds.fa', {})
  #fi
  if "aaURL" in versions[version]:
    files["aa"] = ( versions[version]["aaURL"], 'aa.fa', {})
  if "genomeURLs" in versions[version]:
    for chrID in versions[version]["genomeURLs"]:
      files["chr_%s" % chrID] = (versions[version]["genomeURLs"][chrID], 'chr%s.fa.gz' % chrID, {})
    #efor
  #fi
  return { k : (u, 'genome_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listGenomes():
  print("Available versions:")
  for genome in versions:
    print(" * %s" % genome)
  #efor
#edef

###############################################################################

class Genome(fm.FileManager):

  version = None
  where    = None
  fileIndex = None

  def __init__(self, version, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["gff", "cds", "aa"] + [ ("genome", chrID) for chrID in versions[version]["chr"] ], **kwargs)
    self.version = version

    self.addStrFunction( lambda s: "Genome : %s" % s.version )

    self.gff    = rm.GFF3ResourceManager(self, "gff")
    self.genome = { chrID: rm.FastaResourceManager(self, "chr_%s" % chrID) for chrID in versions[self.version]["chr"] }
    self.cds    = rm.FastaResourceManager(self, "cds")
    self.aa     = rm.FastaResourceManager(self, "aa")
  #edef

  ###############################################################################

#eclass


