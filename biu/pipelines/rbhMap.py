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
  __output = None
  __outputMap = None

  def __init__(self, fasta1, fasta2, config=None, **kwargs):
    Pipeline.__init__(self, snakemakeFile, **kwargs)

    smConfig = self.__defaultConfig
    if isinstance(config, dict):
      smConfig.update(config)
    #fi

    fileName1 = self.__writeTemporaryFile(fasta1)
    fileName2 = self.__writeTemporaryFile(fasta2)

    smConfig["genomes"] = { "A" : { "fasta" : fileName1, "is_prot" : 1 if (fasta1.primaryType == formats.Sequence.PROTTYPE) else 0 },
                            "B" : { "fasta" : fileName2, "is_prot" : 1 if (fasta2.primaryType == formats.Sequence.PROTTYPE) else 0 } }

    self.setConfig(smConfig)
    self.run(["output"])
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

    if self.__outputMap is None:
      self.__outputMap = formats.TSVMap(self.__getOutputFileName(), 0, 1, delimiter='\t', **kwargs)
    #fi
    return self.__outputMap
  #edef

  def getMappingDetails(self):
    if not(self.success):
      return None
    #fi

    if self.__output is None:
      self.__output = pd.read_csv(self.__getOutputFileName(), sep='\t', names=["from", "to", "evalue", "bitscore" ])
    #fi
    return self.__output
  #edef

#eclass
