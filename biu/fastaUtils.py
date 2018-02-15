import gzip

###############################################################################

class Fasta(object):

  entries = None
  fileName = None

  def __init__(self, entries=None, fileName=None):
    self.entries = {}
    if fileName is not None:
      self.entries.update(loadFasta(fileName))
      self.fileName = fileName
    #fi

    if entries is not None:
      self.entries.update(entries)
    #fi
  #edef

  def __str__(self):
    print ("hello")
    dstr  = "Fasta object\n"
    dstr += " Where: %s\n" % self.fileName
    dstr += " Entries: %d\n" % len(self.entries)
    return dstr
  #edef

  def get(self, seqID, start=None, end=None):
    if seqID not in entries:
      return None
    #fi
    seq = self.entries[seqID]
    if start is None:
      start = 0
    #fi
    if end is None:
      end = len(seq)
    #fi
    return seq[start:end]
  #edef
  
  def set(self, seqID, seq):
    self.entries[seqID] = seq
  #edef

  def __getitem__(self, seqID):
    if seqID in self.entries:
      return self.entries[seqID]
    else:
      print("Error, unknown sequence '%s'" % seqID)
      return ""
    #fi
  #edef

  def __setitem__(self, seqID, seq):
    self.entries[seqID] = seq
  #edef

  def update(self, d):
    self.entries.update(d)
  #edef

  def merge(self, other):
    return Fasta({**self.entries, **other.entries})
  #edef

  def write(self, outfile):
    writeFasta(self.entries, outfile)
  #edef

#eclass

def loadFasta(fastaFile):

  F = {'': 0}

  current_seq = ""
  buffer_seq  = ""
  
  with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
    for line in fd:
      line = line.strip()
      if len(line) == 0:
        continue
      #fi
      if line[0] == '>':
        F[current_seq] = buffer_seq
        current_seq = line[1:].split(' ')[0]
        buffer_seq = ""
      else:
        buffer_seq = buffer_seq + line.strip()
      #fi
  #ewith
  F[current_seq] = buffer_seq
  F.pop("", None)
  return F
#edef

###############################################################################

def writeFasta(fasta, outFile, linelength=80):
  with open(outFile, "w") as ofd:
    for  (name, sequence) in fasta:
      ofd.write(">%s\n" % name)
      ofd.write("%s\n" % '\n'.join([sequence[i:i+linelength] for i in range(0, len(sequence), linelength)]))
    #efor
  #ewith
#edef
