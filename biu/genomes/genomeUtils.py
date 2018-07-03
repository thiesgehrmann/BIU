from ..structures import Dataset
from ..config import settings as settings

from .. import formats

class Genome(Dataset):

  version = None
  where    = None
  fileIndex = None

  def __init__(self, version, fileIndex, **kwargs):

    Dataset.__init__(self, fileIndex, **kwargs)

    self._addStrFunction( lambda s: "Genome : %s" % version )
    if "gff" in fileIndex:
      self._registerObject("gff", formats.GFF3, ["gff"], fileIndex["gff"].path)
    #fi

    if "genome" in fileIndex:
      self._registerObject("genome", formats.Fasta, ["genome"], fileIndex["genome"].path, seqType=formats.Sequence.DNATYPE)
    #fi

    if "cds" in fileIndex:
      self._registerObject("cds", formats.Fasta, ["cds"], fileIndex["cds"].path, seqType=formats.Sequence.DNATYPE)
    #fi

    if "aa" in fileIndex:
      self._registerObject("aa", formats.Fasta, ["aa"], fileIndex["aa"].path, seqType=formats.Sequence.PROTTYPE)
    #fi

    if 'ids' in fileIndex:
      self._registerObject('ids', formats.Map, ['ids'], fileIndex["ids"].path, delimiter='\t')
    #fi

    if 'orthology' in fileIndex:
      self._registerObject('orthology', formats.Map, ['orthology'], fileIndex['orthology'].path, delimiter='\t')
    #fi
  #edef

  ###############################################################################

  def seq(self, ID):
    return self.gff.seq(ID, self.genome)
  #edef

#eclass
