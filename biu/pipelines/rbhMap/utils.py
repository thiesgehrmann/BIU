from collections import namedtuple
import csv


class BlastHitType(object):
  BlastHit = None
  def __init__(self):
    self.BlastHit = namedtuple("BlastHit", "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
  def setFields(self, blastFields):
    self.BlastHit = namedtuple("BlastHit", blastFields)
  #edef
  def getFields(self):
    return self.BlastHit._fields

  def blastFieldTypeConversion(self, fieldName, value):
    typeMap = { "qseqid" : lambda x: str(x),
                "sseqid" : lambda x: str(x),
                "pident" : lambda x: float(x),
                "length" : lambda x: int(x),
                "mismatch" : lambda x: int(x),
                "gapopen": lambda x: int(x),
                "qstart": lambda x: int(x),
                "qend": lambda x: int(x),
                "sstart": lambda x: int(x),
                "send": lambda x: int(x),
                "evalue": lambda x: float(x),
                "bitscore": lambda x: float(x),
                "slen": lambda x: int(x),
                "qlen": lambda x: int(x) }
    return typeMap[fieldName](value)
  #edef
#eclass

blastHitType = BlastHitType()

def parseBlastHit(row):

  typedRow = [ blastHitType.blastFieldTypeConversion(field, value) for (field, value) in zip(blastHitType.getFields(), row) ]
  return blastHitType.BlastHit(*typedRow)
#edef

###############################################################################

def blastHit2Row(h):
  return [str(x) for x in [h[i] for i in range(len(h._fields))]]
#edef

###############################################################################

def readBlastFile(filename, fields="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
  blastHitType.setFields(fields)
  hits    = []
  # Read the BLAST hits
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t", quotechar="\"");
    for row in reader:
      hits.append(parseBlastHit(row))
    #efor
  #ewith
  return hits
#edef

###############################################################################

def writeBlastFile(H, filename):
  with open(filename, "w") as ofd_hits:
    for h in H:
      ofd_hits.write('\t'.join(blastHit2Row(h)) + '\n')
    #efor
  #ewith
#edef

###############################################################################

def indexListBy(L, key=lambda x: x[0]):
  G = {}
  for item in L:
    k = key(item)
    if k not in G:
      G[k] = []
    #fi
    G[k].append(item)
  #efor
  return G
#edef
