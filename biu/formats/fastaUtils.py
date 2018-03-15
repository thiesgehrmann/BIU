import gzip

from .. import utils

from .seqUtils import Sequence

###############################################################################

class Fasta(object):

  entries = None
  fileName = None

  def __init__(self, data, entries=None, fileName=None, seqType=Sequence.DNATYPE):
    if isinstance(data, str):
      utils.dbm("Fasta input source is file")
      self.entries = Fasta.loadFasta(data, seqType)
      self.fileName = data
    else:
      if isinstance(data, dict):
        utils.dbm("Fasta input source is dictionary of sequences.")
        self.entries = data
      elif hasattr(data, "__iter__"):
        utils.dbm("Fasta input source is a list of sequences.")
        self.entries = { s.name : s if isinstance(s, Sequence) else Sequence(str(i), s, seqType)  for i,s in enumerate(data) }
      #fi
    #fi
  #edef

  def __str__(self):
    dstr  = "Fasta object\n"
    dstr += " Where: %s\n" % self.fileName
    dstr += " Entries: %d\n" % len(self.entries)
    dstr += " Primary type: %s\n" % self.primaryType
    return dstr
  #edef

  @property
  def primaryType(self):
    for k in self.entries:
      return self.entries[k].seqType
    #efor
    return Sequence.UNKNOWNTYPE
  #edef

  def keys(self):
    return self.entries.keys()
  #edef

  def items(self):
    return self.entries.items()
  #edef

  def values(self):
    return self.entries.values()
  #edef

  def __iter__(self):
    self._iterKeys = list(self.keys())
    return self
  #edef

  def __next__(self):
    if len(self._iterKeys) == 0:
      raise StopIteration
    #fi
    v = self._iterKeys.pop()
    return v
  #edef

  def __getitem__(self, seqID):
    if seqID in self.entries:
      return self.entries[seqID]
    else:
      utils.error("Unknown sequence '%s'" % seqID)
      return None
    #fi
  #edef

  def __setitem__(self, seqID, seq):
    if not(isinstance(seq, Sequence)):
      seq = Sequence(seqID, seq, self.primaryType)
    #fi
    self.entries[seqID] = seq
  #edef

  def update(self, d):
    self.entries.update(d)
  #edef

  def merge(self, other):
    return Fasta({**self.entries, **other.entries})
  #edef

  def write(self, outfile):
    Fasta.writeFasta(self.entries, outfile)
  #edef

  ###############################################################################

  @staticmethod
  def loadFasta(fastaFile, seqType=Sequence.DNATYPE):
  
    F = {'': 0}
  
    current_seq = ""
    current_seq_full = ""
    buffer_seq  = ""
    
    with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
      for line in fd:
        line = line.strip()
        if len(line) == 0:
          continue
        #fi
        if line[0] == '>':
          F[current_seq] = Sequence(current_seq, buffer_seq, seqType, current_seq_full)
          current_seq = line[1:].split(' ')[0]
          current_seq_full = line[1:]
          buffer_seq = ""
        else:
          buffer_seq = buffer_seq + line.strip()
        #fi
    #ewith
    F[current_seq] = Sequence(current_seq, buffer_seq, seqType, current_seq_full)
    F.pop("", None)
    return F
  #edef
  
  ###############################################################################

  @staticmethod
  def writeFasta(fasta, outFile, linelength=80):
    with open(outFile, "w") as ofd:
      for  seqid in fasta:
        seq = fasta[seqid]
        name = seq.fullName
        seq  = seq.seq
        ofd.write(">%s\n" % name)
        ofd.write("%s\n" % '\n'.join([seq[i:i+linelength] for i in range(0, len(seq), linelength)]))
      #efor
    #ewith
  #edef

#eclass
