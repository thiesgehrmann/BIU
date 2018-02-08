import csv
from collections import namedtuple
import tabix

from . import fileManager as fm

###############################################################################
versions = { "GRCh37" : {
  "caddURL"      : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz",
  "caddTabixURL" : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi"
  }
}

def listVersions():
  print("Available versions:")
  for v in versions:
    print(" * %s" % v)
#edef

###############################################################################

class CADD(fm.FileManager):

  version = None

  def __init__(self, version=list(versions.keys())[0], where='./'):
    fm.FileManager.__init__(self, where + '/%s' % version)
    self.version = version
    self.fileIndex = self.__urlFileIndex()
    self.str_functions.append(lambda s: "Version: %s" % self.version)

    def loadedObjects(s):
      dstr = 'Objects:\n'
      dstr += " * [%s] _source\n" % ('X' if s._source is not None else ' ')
      return dstr
    #edef

    self.str_functions.append(lambda s: loadedObjects(s))
  #edef

  def __urlFileIndex(self):
    files = {}

    files["tsv"] = (versions[self.version]["caddURL"], self.where + '/cadd.tsv.bgz', {})
    files["tsv_tbi"] = (versions[self.version]["caddTabixURL"], self.where + '/cadd.tsv.bgz.tbi', {})

    return files
  #edef

  #############################################################################

  caddFields = [ "chrom", "pos", "ref", "alt", "rawscore", "phred" ]
  caddEntry = namedtuple("CADDEntry", caddFields);
  
  _source = None

  def _requireSource(self):
    if self._source is None:
      if not(self.satisfyRequiredFiles(["tsv"])):
        return False
      #fi
      self._source = tabix.open(self.getFileName("tsv"))
    #fi
    return True
  #fi

  def query(self, chromosome, start, end=None, alt=None):
    if not(self._requireSource()):
      return None
    #fi

    qres = utils.tabixQueryWrapper(self._source, chromosome, start, (start if (end is None) else end))
    qres = [ self.caddEntry(*r) for r in qres ]
    if (alt is None) and (end is None):
      resPhred = {}
      for res in qres:
        resPhred[res.alt] = res.phred
      #efor
      return resPhred
    #fi
    if (alt is None) and (end is not None):
      resPhred = {}
      for res in qres:
        resPhred[(int(res.pos),res.alt)] = res.phred
      #efor
      return resPhred
    #fi
    if (alt is not None) and (end is not None):
      resPhred = {}
      for res in [r for r in qres if r.alt == alt]:
        resPhred[int(res.pos)] = res.phred
      #efor
      return resPhred
    #fi
    if (end is None) and (alt is not None):
      relRes = [ r for r in qres if r.alt == alt ]
      if len(relRes) != 1:
        return None
      else:
        return relRes[0].phred
      #fi
    #fi
  #edef

#eclass
    
