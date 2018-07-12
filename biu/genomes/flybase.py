from .genomeUtils import Genome

from .. import utils
from .. import formats
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

  def __init__(self, release='FB2018_03', organism="dmel_r6.22", where=None, **kwargs):
    version = "flybase_%s.%s" % (release, organism)
    fileIndex = self.__genFileIndex(version, release, organism, where=where)
    Genome.__init__(self, version, fileIndex, **kwargs)
  #edef

  def __genFileIndex(self, version, release, organism, where):
    files = {}
    finalPath = '%s/genomes/flybase/%s' % ( (settings.getDataDir() if where is None else where), version)

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

    def genIDS():
      def idMapFunc(inFile, outFile):
        fasta = formats.Fasta(inFile)
        with open(outFile, 'w') as ofd:
          ofd.write('\t'.join([ 'gene', 'transcript', 'protein', 'symbol' ]) + '\n')
          for seq in fasta:
            data = dict([ (s[0], s[1]) for s in  [ p.split('=') for p in fasta[seq].fullName.split(' ') if '=' in p ] if s[0] in ['name', 'parent'] ])
            ofd.write('\t'.join([data.get('parent', ',').split(',')[0], data.get('parent', ',').split(',')[1].replace(';', ''), seq, data.get('name', '') ]) + '\n')
          #efor
        #ewith
        return 0
      #edef
      uri = [ line for line in conn.nlst("/releases/%s/%s/fasta" % (release, organism)) if 'all-translation' in line ]
      if len(uri) > 0:
        return utils.Acquire(where=where).curl("ftp://ftp.flybase.net/%s" % uri[0]).gunzip().func(idMapFunc).finalize('%s/ids.tsv' % finalPath)
      #fi
      return None
    #edef

  

    files["gff"]    = genGFF3()
    files["genome"] = genGenome()
    files["cds"]    = genCDS()
    files["aa"]     = genAA()
    files['ids']    = genIDS()

    for f in files:
      if f is None:
        del files[f]
      #fi
    #efor

    conn.close()
    return files
  #edef
#eclass
