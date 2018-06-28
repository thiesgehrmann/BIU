from .genomeUtils import Genome

from .. import utils
from ..config import settings

import csv
import re
from collections import namedtuple

class Wormbase(Genome):

  @classmethod
  def organisms(cls):
    print("Organisms in Wormbase")
    for organism in cls.__organisms():
      print(" * %s" % organism)
    #efor
  #edef

  @staticmethod
  def __organisms():
    speciesTuple = namedtuple('wormbaseSpecies', ['species', 'bioproject', 'ftp', 'genome', 'gff3', 'aa'])
    
    ao = utils.Acquire().curl('https://wormbase.org/rest/widget/index/all/all/downloads')
    
    pat = re.compile("<a\s+href=\"(ftp.*?)\"[^>]*>.*?</a>",re.M|re.DOTALL|re.I)
    htmlString = pat.sub('\\1', ''.join(open(ao.acquire(), 'r').readlines()))
    rows = [ [ f.strip().replace(' ', '_') for f in row ] for row in utils.html.table2csv(''.join(htmlString)) ]
    rows = [ speciesTuple(*row) for row in rows if len(row) == 6 ]
    return { s.species : s for s in rows }
  #edef

  def __init__(self, organism="Caenorhabditis_elegans", where=None):
    version = "wormbase_%s" % (organism)
    fileIndex = self.__genFileIndex(organism, where=where)
    Genome.__init__(self, version, fileIndex)
  #edef

  def __genFileIndex(self, organism, where):
    organisms = self.__organisms()
    if organism not in organisms:
      utils.msg.warning("Organism '%s' not in Wormbase" % (organism))
      return {}
    #fi
    species = organisms[organism]
    files = {}
    finalPath = '%s/wormbase_%s' % ( (settings.getWhere() if where is None else where), organism)

    files["gff"]    = utils.Acquire(where=where).curl(species.gff3).gunzip().finalize('%s/genes.gff3' % finalPath)
    files["genome"] = utils.Acquire(where=where).curl(species.genome).gunzip().finalize('%s/genome.fasta' % finalPath)
    files["aa"]     = utils.Acquire(where=where).curl(species.aa).gunzip().finalize('%s/aa.fasta' % finalPath)

    return files
  #edef

#eclass
    
