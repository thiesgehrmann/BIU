from .genomeUtils import Genome

from .. import utils
from ..config import settings


import ftplib

class Ensembl(Genome):

  @staticmethod
  def organisms(release=92, grch37=False):
    basedir = '/pub/grch37' if grch37 else '/pub'
    conn      = ftplib.FTP("ftp.ensembl.org")
    conn.login()
    organisms = [ line.split('/')[-1] for line in conn.nlst("%s/release-%d/gff3" % (basedir, release)) ]
    print("Organisms in Ensembl, release %d:" % release)
    for organism in organisms:
      print(" * %s" % organism)
    #efor
  #edef

  def __init__(self, release=92, organism="homo_sapiens", grch37=False, where=None):
    self.__grch37 = grch37
    basedir = '/pub/grch37' if grch37 else '/pub'
    version = "ensembl_%s%s.%s" % ('grch37.' if self.__grch37 else '', str(release), organism)
    fileIndex = self.__genFileIndex(version, basedir, release, organism, where=where)
    Genome.__init__(self, version, fileIndex)
  #edef

  def __genFileIndex(self, version, basedir, release, organism, where):
    files = {}
    finalPath = '%s/%s' % ( (settings.getWhere() if where is None else where), version)
    conn      = ftplib.FTP("ftp.ensembl.org")
    conn.login()
    organisms = [ line.split('/')[-1] for line in conn.nlst("%s/release-%d/gff3" % (basedir, release)) ]

    if organism not in organisms:
      utils.msg.warning("Organism '%s' not in release %d" % (organism, release))
      return {}
    #fi

    def genGFF3():
      # This is ugly, but the most stable way I was able to find the correct GFF3 file.
      # in GRCH37, it didnt match the obvious pattern...
      uri = [ u[0] for u in  sorted([ (line, len(line.split('.'))) for line in conn.nlst("%s/release-%d/gff3/%s" % (basedir, release, organism)) if (len(line.split('.')) > 3) ],
                                    key=lambda x:x[1]) ]
      #uri = [ line for line in conn.nlst("%s/release-%d/gff3/%s" % (basedir, release, organism)) if '%d.gff3.gz' % release in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.ensembl.org/%s" % uri[0]).gunzip().finalize('%s/genes.gff3' % finalPath)
      #fi
      return None
    #edef

    def genGenome():
      uri = [ line for line in conn.nlst("%s/release-%d/fasta/%s/dna" % (basedir, release, organism)) if 'dna.chromosome' in line ]
      aos = [ utils.Acquire(where=where).curl("ftp://ftp.ensembl.org/%s" % f) for f in uri ]
      return utils.Acquire(where=where).merge(aos, method='zcat').finalize('%s/dna.fasta' % finalPath)
    #edef

    def genCDS():
      uri = [ line for line in conn.nlst("%s/release-%d/fasta/%s/cds" % (basedir, release, organism)) if 'fa.gz' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.ensembl.org/%s" % uri[0]).gunzip().finalize('%s/cds.fa' % finalPath)
      #fi
      return None
    #edef

    def genAA():
      uri = [ line for line in conn.nlst("%s/release-%d/fasta/%s/pep" % (basedir, release, organism)) if 'all.fa.gz' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.ensembl.org/%s" % uri[0]).gunzip().finalize('%s/aa.fa' % finalPath)
      #fi
      return None
    #edef

    files["gff"]    = genGFF3()
    files["genome"] = genGenome()
    files["cds"]    = genCDS()
    files["aa"]     = genAA()

    for f in files:
      if f is None:
        del files[f]
      #fi
    #efor

    conn.close()
    return files
  #edef
#eclass
