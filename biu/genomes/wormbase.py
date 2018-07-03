from .genomeUtils import Genome
from .. import formats
from .. import utils
from ..config import settings

import csv
import re
from collections import namedtuple
import ftplib

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
    rows = [ row for row in rows if len(row) == 6 ]

    # Annoyingly, some species have multiple lines, with the same species identifier. Resolve that here.
    renamedRows = {}
    for row in rows:
      origSpecies = row[0]
      species = row[0]
      speciesDupIdx = 0
      while species in renamedRows:
        speciesDupIdx += 1
        species = '%s.%d' % (origSpecies, speciesDupIdx)
      #ewhile
      renamedRows[species] = speciesTuple(species, *row[1:])
    #efor

    return renamedRows
  #edef

  def __init__(self, organism="Caenorhabditis_elegans", where=None, **kwargs):
    version = "wormbase_%s" % (organism)
    fileIndex = self.__genFileIndex(organism, where=where)
    Genome.__init__(self, version, fileIndex, **kwargs)
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
    files['cds']    = utils.Acquire(where=where).curl(species.aa.replace('protein.fa.gz', 'CDS_transcripts.fa.gz')).gunzip().finalize('%s/cds.fasta' % finalPath)

    def idMapFunc(inFile, outFile):
      fasta = formats.Fasta(inFile)
      with open(outFile, 'w') as ofd:
        ofd.write('\t'.join([ 'gene', 'protein', 'peptide', 'uniprot', 'insdc' ]) + '\n')
        for seq in fasta:
          data = dict([ (s[0], s[1]) for s in  [ p.split('=') for p in fasta[seq].fullName.split(' ') if '=' in p ] if s[0] in ['wormpep', 'gene', 'uniprot', 'insdc'] ])
          ofd.write('\t'.join([data.get('gene', ''), seq, data.get('wormpep', ''), data.get('uniprot', ''), data.get('insdc', '') ]) + '\n')
        #efor
      #ewith
      return 0
    #edef
    files['ids'] = utils.Acquire(where=where).curl(species.aa).gunzip().func(idMapFunc).finalize('%s/map.tsv' % finalPath)

    def orthologyMapFunc(inFile, outFile):
      species = set([])
      mapping = {}
      gene = None
      with open(inFile, 'r') as ifd:
        while True:
          line = ifd.readline()
          if line == '':
            break
          #fi
          line = line.strip()
          if line[0] == '#':
            continue
          elif line[0] == '=':
            geneLine = ifd.readline().strip()
            gene = geneLine.split('\t')[0]
            mapping[gene] = {}
          else:
            lineSpecies, lineGene, lineSource = line.split('\t')[:3]
            lineSpecies = lineSpecies.replace(' ', '_')
            species.add(lineSpecies)
            mapping[gene][lineSpecies] = lineGene
          #fi
        #efor
      #ewith
      species = sorted(list(species))
      with open(outFile, 'w') as ofd:
        ofd.write('gene\t%s\n' % ('\t'.join(species)))
        for gene in mapping:
          ofd.write('%s\t%s\n' % (gene, '\t'.join([ mapping[gene].get(s, '') for s in species ])))
        #efor
      #ewith
      return 0
    #edef

    conn = ftplib.FTP("ftp.wormbase.org")
    conn.login()
    uri = [ line for line in conn.nlst('/%s/annotation' % '/'.join(species.aa.split('/')[3:-1])) if 'orthologs' in line ]
    if len(uri) > 0:
      files['orthology'] = utils.Acquire(where=where).curl('ftp://ftp.wormbase.org/%s' % uri[0]).gunzip().func(orthologyMapFunc).finalize('%s/orthology.tsv' % finalPath)
    #fi

    return files
  #edef

#eclass
    
