from .. import utils

###############################################################################

bases = ['t', 'c', 'a', 'g'];
codons = [a+b+c for a in bases for b in bases for c in bases];
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
codon_table = dict(zip(codons, amino_acids));
codon_table_rev = dict(zip(amino_acids, codons))
revcompDict = { 'A' : 'T', 'T' : 'A', 'a' : 't', 't' : 'a',
                'C' : 'G', 'G' : 'C', 'c' : 'g', 'g' : 'c',
                'R' : 'Y', 'Y' : 'R', 'r' : 'y', 'y' : 'r',
                'K' : 'M', 'M' : 'K', 'k' : 'm', 'm' : 'k',
                'S' : 'W', 'W' : 'S', 's' : 'w', 'w' : 's',
                'B' : 'V', 'V' : 'B', 'b' : 'v', 'v' : 'b',
                'D' : 'H', 'H' : 'D', 'd' : 'h', 'h' : 'd',
                'N' : 'N', 'X' : 'X', 'n' : 'n', 'x' : 'x' }

###############################################################################

class Sequence(object):

  DNATYPE = 'dna'
  PROTTYPE = 'prot'
  UNKNOWNTYPE = 'unknown'

  __name = None
  __fullName = None
  __seq = None
  __seqType = None

  def __init__(self, name, seq, seqType = DNATYPE, fullName=None):
    self.__name = name
    self.__fullName = fullName
    self.__seq = seq
    self.__seqType = seqType
  #edef

  def __str__(self):
    return self.__seq
  #edef

  def __getitem__(self, s):
    return self.seq[s]
  #edef

  @property
  def name(self):
    return self.__name
  #edef

  @property
  def fullName(self):
    if self.__fullName is None:
      return self.name
    else:
      return self.__fullName
    #fi
  #edef

  @property
  def seq(self):
    return self.__seq
  #edef

  @property
  def seqType(self):
    return self.__seqType
  #edef

  ###############################################################################
  
  def translate(self):
    if len(self.__seq) % 3 != 0:
      utils.warning("The sequence you are trying to translate is not divisible by 3.")
    #fi
    aa = '';
    for i in range(0, len(self.__seq), 3):
      cod = self.__seq[i:i+3].lower()
      aa += codon_table[cod] if (cod in codon_table) else '*'
    #efor
    return Sequence(self.__name, aa, self.PROTTYPE, self.__fullName);
  #edeif

  def reverseTranslate(self):
    if self.__seqType == self.PROTTYPE:
      nt = ''.join([ codon_table_rev[a] for a in self.__seq])
      return Sequence(self.__name, nt, self.DNATYPE, self.__fullName)
    else:
      utils.error("Cannot do reverse translation on DNA type.")
      return None
    #fi
  #edef
  
  ###############################################################################
  
  def revcomp(self):
    if self.__seqType == self.DNATYPE:
      return Sequence(self.name, ''.join([ revcompDict[b] for b in self.__seq[::-1] ]), self.seqType, self.__fullName)
    else:
      utils.error("Cannot do reverse complement on PROT type.")
      return None
    #fi
  #edef
  
  ###############################################################################

#eclass
