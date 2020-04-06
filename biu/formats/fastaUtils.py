import gzip

from .. import utils

from .seqUtils import Sequence

###############################################################################

class Fasta(object):

  __slots__ = [ '__entries', '__fileName', '__iterKeys' ]

  #__entries = None
  #__fileName = None
  #__iterKeys = None

  def __init__(self, data, seqType=Sequence.DNATYPE):
    self.__fileName = None

    if isinstance(data, str):
      utils.dbm("Fasta input source is file")
      self.__entries = Fasta.loadFasta(data, seqType)
      self.__fileName = data
    else:
      if isinstance(data, dict):
        utils.dbm("Fasta input source is dictionary of sequences.")
        self.__entries = data
      elif hasattr(data, "__iter__"):
        utils.dbm("Fasta input source is a list of sequences.")
        self.__entries = {}
        for i,s in enumerate(data):
            if isinstance(s, Sequence):
                self.__entries[s.name] = s
            else:
                self.__entries[str(i)] = Sequence(str(i), s, seqType)
            #fi
        #efor
      #fi
    #fi
  #edef

  def __str__(self):
    dstr  = "Fasta object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    dstr += " Entries: %d\n" % len(self.__entries)
    dstr += " Primary type: %s\n" % self.primaryType
    return dstr
  #edef
    
  def __repr__(self):
    return str(self)
  #edef

  @property
  def primaryType(self):
    for k in self.__entries:
      return self.__entries[k].seqType
    #efor
    return Sequence.UNKNOWNTYPE
  #edef

  def keys(self):
    return self.__entries.keys()
  #edef

  def items(self):
    return self.__entries.items()
  #edef

  def values(self):
    return self.__entries.values()
  #edef

  def __contains(self, k):
    return (k in self.__entries)
  #edef

  def __iter__(self):
    self.__iterKeys = list(self.__entries.keys())
    return self
  #edef

  def __next__(self):
    if len(self.__iterKeys) == 0:
      raise StopIteration
    #fi
    v = self.__iterKeys.pop()
    return v
  #edef

  def __getitem__(self, seqID):
    if seqID in self.__entries:
      return self.__entries[seqID]
    else:
      utils.error("Unknown sequence '%s'" % seqID)
      return None
    #fi
  #edef

  def __setitem__(self, seqID, seq):
    if not(isinstance(seq, Sequence)):
      seq = Sequence(seqID, seq, self.primaryType)
    #fi
    self.__entries[seqID] = seq
  #edef

  def update(self, d):
    self.__entries.update(d)
  #edef

  def merge(self, other):
    return Fasta({**self.__entries, **other.__entries})
  #edef

  def write(self, outfile):
    Fasta.writeFasta(self.__entries, outfile)
  #edef

  ###############################################################################

  @staticmethod
  def loadFasta(fastaFile, seqType=Sequence.DNATYPE):
  
    F = {'': 0}
  
    current_seq = ""
    current_seq_full = ""
    buffer_seq  = ""
   
    with utils.gzopen(fastaFile, 'rt', encoding='UTF-8') as fd: 
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
          try:
            buffer_seq = buffer_seq + line
          except TypeError:
            print(type(buffer_seq))
            
            print(type(line))
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
