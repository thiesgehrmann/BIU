from .genomeUtils import Genome

from .. import utils
from ..config import settings

import ftplib

class Flybase(Genome):

  @classmethod
  def organisms(cls, release='FB2018_03'):
    print("Organisms in Flybase, release %s:" % release)
    for organism in cls.__organisms(release):
      print(" * %s" % organism)
    #efor
  #edef

  @staticmethod
  def __organisms(release):
    conn = ftplib.FTP("ftp.flybase.net")
    conn.login()
    organisms = [ line.split('/')[-1] for line in conn.nlst("releases/%s" % release) if '_r' in line ]
    return organisms
  #edef

  def __init__(self, release='FB2018_03', organism="dmel_r6.22", where=None):
    version = "flybase_%s.%s" % (release, organism)
    fileIndex = self.__genFileIndex(version, release, organism, where=where)
    Genome.__init__(self, version, fileIndex)
  #edef

  def __genFileIndex(self, version, release, organism, where):
    files = {}
    finalPath = '%s/%s' % ( (settings.getWhere() if where is None else where), version)

    conn = ftplib.FTP("ftp.flybase.net")
    conn.login()

    if organism not in self.__organisms(release):
      utils.msg.warning("Organism '%s' not in release %d" % (organism, release))
      return {}
    #fi

    def genGFF3():
      uri = [ line for line in conn.nlst("/releases/%s/%s/gff" % (release, organism)) if 'all-filtered' in line ]
      #uri = [ line for line in conn.nlst("%s/release-%d/gff3/%s" % (basedir, release, organism)) if '%d.gff3.gz' % release in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.flybase.net/%s" % uri[0]).gunzip().finalize('%s/genes.gff3' % finalPath)
      #fi
      return None
    #edef

    def genGenome():
      uri = [ line for line in conn.nlst("/releases/%s/%s/fasta" % (release, organism)) if 'all-chromosome' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.flybase.net/%s" % uri[0]).finalize('%s/dna.fasta' % finalPath)
      #fi
      return None
    #edef

    def genCDS():
      uri = [ line for line in conn.nlst("/releases/%s/%s/fasta" % (release, organism)) if 'all-transcript' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.flybase.net/%s" % uri[0]).gunzip().finalize('%s/cds.fa' % finalPath)
      #fi
      return None
    #edef

    def genAA():
      uri = [ line for line in conn.nlst("/releases/%s/%s/fasta" % (release, organism)) if 'all-translation' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.flybase.net/%s" % uri[0]).gunzip().finalize('%s/aa.fa' % finalPath)
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
