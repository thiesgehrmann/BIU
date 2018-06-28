from .pipeline import Pipeline
from .. import formats

import inspect, os
import pandas as pd

snakemakeFile = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/rbhMap/Snakefile'

class RBHMap(Pipeline):

  __defaultConfig = {
    "e_threshold" : 10e-10,
    "blast_fields" : "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen",
    "output_file_name" : "mapping.tsv"
  }
  __output = {}

  def __init__(self, fasta1, fasta2, config={}, **kwargs):
    Pipeline.__init__(self, snakemakeFile, {**self.__defaultConfig, **config}, **kwargs)

    fileName1 = self.__writeTemporaryFile(fasta1)
    fileName2 = self.__writeTemporaryFile(fasta2)

    genomes = { "A" : { "fasta" : fileName1, "is_prot" : 1 if (fasta1.primaryType == formats.Sequence.PROTTYPE) else 0 },
                "B" : { "fasta" : fileName2, "is_prot" : 1 if (fasta2.primaryType == formats.Sequence.PROTTYPE) else 0 } }

    self.setConfig(genomes=genomes)
    if self.autorun:
      self.run(["output"])
    #fi
  #edef

  def __writeTemporaryFile(self, fasta, hashName=True):
    fileName, exists = self._generateInputFileName([ fasta[s].seq for s in fasta ])

    if not(exists):
      fasta.write(fileName)
    #fi

    return fileName
  #edef 

  def __getOutputFileName(self):
    return '%s/%s' % (self.config["outdir"], self.config["output_file_name"]) 
  #edef

  def getMapping(self, **kwargs):
    if not(self.success):
      return None
    #fi

    if "map" not in self.__output:
      self.__output["map"] = formats.TSVMap(self.__getOutputFileName(), 0, 1, delimiter='\t', **kwargs)
    #fi
    return self.__output["map"]
  #edef

  def getMappingDetails(self):
    if not(self.success):
      return None
    #fi

    if "tbl" not in self.__output:
      self.__output["tbl"] = pd.read_csv(self.__getOutputFileName(), sep='\t', names=["from", "to", "evalue", "bitscore" ])
    #fi
    return self.__output["tbl"]
  #edef

#eclass
