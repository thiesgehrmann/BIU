from ..structures import Pipeline
from .. import formats
from .. import utils

import inspect, os

pd = utils.py.loadExternalModule("pandas")

snakemakeFile = { False: os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/blast/Snakefile',
                  True: os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/blast/diamond.Snakefile'}

class Blast(Pipeline):

  __defaultConfig = {
    "e_threshold" : 10e-10,
    "blast_fields" : "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen",
    "max_target_seqs" : 10,
    "output_file_name" : "alignment.tsv"
  }
  __output = None

  def __init__(self, fasta1, fasta2, config={}, diamond=False, **kwargs):
    Pipeline.__init__(self, snakemakeFile[diamond], {**self.__defaultConfig, **config}, **kwargs)

    fileName1 = fileName2 = None
    if diamond:
      fileName1 = self.__writeTemporaryFile(formats.Fasta([ fasta1[s].translate() if (fasta1[s].seqType == formats.Sequence.DNATYPE ) else fasta1[s] for s in fasta1 ]))
      fileName2 = self.__writeTemporaryFile(formats.Fasta([ fasta2[s].translate() if (fasta2[s].seqType == formats.Sequence.DNATYPE ) else fasta2[s] for s in fasta2 ]))
    else:
      fileName1 = self.__writeTemporaryFile(fasta1)
      fileName2 = self.__writeTemporaryFile(fasta2)
    #fi

    genomes = { "A" : { "fasta" : fileName1, "is_prot" : 1 if (fasta1.primaryType == formats.Sequence.PROTTYPE) else 0 },
                "B" : { "fasta" : fileName2, "is_prot" : 1 if (fasta2.primaryType == formats.Sequence.PROTTYPE) else 0 } }

    self.setConfig(genomes=genomes, diamond=(1 if diamond else 0))
    if self.autorun:
      self.run(["output"])
    #fi
    self.__output = {}
  #edef

  def __writeTemporaryFile(self, fasta):
    fileName, exists = self._generateInputFileName([ fasta[s].seq for s in fasta ])

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

    if "map" not in self.__output:
      self.__output["map"] = formats.Map(self.__getOutputFileName(), names=self.config["blast_fields"].split(' '), onlyIndex=['qseqid', 'sseqid'], delimiter='\t', **kwargs)
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
