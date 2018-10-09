from ..structures import Pipeline
from .. import formats
from .. import utils

import inspect, os

pd = utils.py.loadExternalModule("pandas")

snakemakeFile = { "blast": os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/blast/Snakefile',
                  "diamond": os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/blast/diamond.Snakefile',
                  "blat" : os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/blast/blat.Snakefile' }

class Blast(Pipeline):
  """
  Align using Blast, Diamond, or BLAT.
  
  """

  __defaultConfig = {
    "e_threshold" : 10e-10,
    "blast_fields" : "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen",
    "options" : "",
    "output_file_name" : "alignment.tsv"
  }
  __output = None

  def __init__(self, database, query, config={}, diamond=False, blat=False, **kwargs):
    config['algorithm'] = algorithm = 'blat' if blat else ('diamond' if diamond else 'blast')
    Pipeline.__init__(self, snakemakeFile[config['algorithm']], {**self.__defaultConfig, **config}, **kwargs)

    fileName1 = fileName2 = None
    if diamond:
      fileName1 = self.__writeTemporaryFile(formats.Fasta([ database[s].translate() if (database[s].seqType == formats.Sequence.DNATYPE ) else database[s] for s in database ]))
      fileName2 = self.__writeTemporaryFile(formats.Fasta([ query[s].translate() if (query[s].seqType == formats.Sequence.DNATYPE ) else query[s] for s in query ]))
    else:
      fileName1 = self.__writeTemporaryFile(database)
      fileName2 = self.__writeTemporaryFile(query)
    #fi

    genomes = { "A" : { "fasta" : fileName1, "is_prot" : 1 if (database.primaryType == formats.Sequence.PROTTYPE) else 0 },
                "B" : { "fasta" : fileName2, "is_prot" : 1 if (query.primaryType == formats.Sequence.PROTTYPE) else 0 } }

    self.setConfig(genomes=genomes, diamond=(1 if diamond else 0))
    if self.autorun:
      self.run(["output"])
    #fi
    self.__output = {}
  #edef

  def __writeTemporaryFile(self, fasta):
    fileName, exists = self._generateInputFileName([ '>%s\n%s' % (s, fasta[s].seq) for s in fasta ])

    if not(exists):
      fasta.write(fileName)
    #fi

    return fileName
  #edef 

  def __getOutputFileName(self):
    print((self.config["outdir"], self.config["output_file_name"]))
    return '%s/%s' % (self.config["outdir"], self.config["output_file_name"]) 
  #edef

  def getResult(self, **kwargs):
    if not(self.success):
      return None
    #fi

    if "tbl" not in self.__output:
      if self.config['algorithm'] == 'blat':
        names = [ 'match', 'mismatch', 'repmatch', 'ns',
                  'gapopen', 'q_gap_bases', 
                  't_gap_count', 't_gap_bases', 'strand',
                  'qseqid', 'qlen', 'qstart', 'qend',
                  'sseqid', 'slen', 'sstart', 'send',
                  'block_count', 'block_sizes',
                  'q_starts', 't_starts' ]
        tbl = pd.read_csv(self.__getOutputFileName(), sep='\t', names=names, skiprows=5)
        if tbl.shape[1] > 0:
          tbl["score"], tbl["pident"] = list(zip(*tbl.apply(lambda r: blatScores(tuple(r)), axis=1)))
          tbl["length"] = tbl["match"] + tbl["mismatch"] + tbl["q_gap_bases"]
        else:
          df['score'] = None
          df['pident'] = None
          df['length'] = None
        #fi
        tbl = tbl[["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "score", "slen", "qlen"]]
        self.__output["tbl"] = tbl
      else:
        self.__output["tbl"] = pd.read_csv(self.__getOutputFileName(), sep='\t', names=self.config["blast_fields"].split(' '))
    #fi
    return self.__output["tbl"]
  #edef

#eclass

def blatScores(psl_row):
  """
   Takes as input the BLAT PSL output rows, and produces the necessary scores
   Ported from perl at: https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
  """
  
  import math
  def pslCalcMilliBad(sizeMul, qEnd, qStart, tEnd, tStart, qNumInsert, tNumInsert,
                      matches, repMatches, misMatches, isMrna):
      """
      Calculate a BLAT score.
      
      """
      milliBad = 0
      qAliSize = sizeMul * (qEnd - qStart)
      tAliSize = tEnd - tStart;
      aliSize  = min(tAliSize, qAliSize)
  
      if aliSize <= 0:
          return 0
      #fi
  
      sizeDif = qAliSize - tAliSize
      if sizeDif < 0:
          sizeDif = 0 if isMrna else -sizeDif
      #fi
      insertFactor = qNumInsert
      if not isMrna:
          insertFactor += tNumInsert
      #fi
      total = sizeMul * (matches + repMatches + misMatches)
      if (total != 0):
          roundAwayFromZero = 3 * math.log(1 + sizeDif)
          if roundAwayFromZero < 0:
              roundAwayFromZero = int(roundAwayFromZero - 0.5)
          else:
              roundAwayFromZero = int(roundAwayFromZero + 0.5)
          #fi
          milliBad = (1000 * (misMatches*sizeMul + insertFactor + roundAwayFromZero)) / total
      #fi
      return milliBad
  #edef
  
  def pslIsProtein(strand, tStart, tEnd, tSize, tStarts, blockSizes):
      starts = [ f for f in tStarts.split(',') if f != "" ]
      sizes  = [ f for f in blockSizes.split(',') if f != "" ]
      answer = False
      if strand == '+':
          test = starts[-1] + (3 * sizes[-1])
          answer = (tEnd == test)
      elif strand == '-':
          test = tSize - (starts[-1] + (3 * sizes[-1]))
          answer = (tStart == test)
      #fi
      return answer
  #edef
  
  """
  Takes as input the BLAT PSL output rows, and produces the necessary scores
  """
  (matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert,
   tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd,
   blockCount, blockSizes, qStarts, tStarts) = psl_row
  sizeMul = 3 if pslIsProtein(strand, tStart, tEnd, tSize, tStarts, blockSizes) else 1
  pslScore = sizeMul * (matches + ( repMatches >> 1) ) - sizeMul * misMatches - qNumInsert - tNumInsert;
  milliBad = int(pslCalcMilliBad(sizeMul, qEnd, qStart, tEnd, tStart, 
                                 qNumInsert, tNumInsert, matches, repMatches, misMatches, 1))
  percentIdentity = 100.0 - milliBad * 0.1;
  return pslScore, percentIdentity

#edef
