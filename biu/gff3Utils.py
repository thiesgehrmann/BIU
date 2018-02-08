import errno
import os
import csv
import gzip
from collections import namedtuple


import gzip

###############################################################################

gff3Entry = namedtuple("gff3Entry", "seqid, source, feature, start, end, score, strand, phase, attr")

def parseGFF3entry(fields):

  def attrsplit(attr):
    spl = attr.split('=')
    if len(spl) == 1:
      return (attr, None)
    elif len(spl) > 2:
      return (spl[0], '='.join(spl[1:]))
    else:
      return (spl[0], spl[1])
    #fi
  #edef

  (seqid, source, feature, start, end, score, strand, phase, attr) = fields
  attr = dict([attrsplit(x.strip()) for x in attr.split(";") ])
  return gff3Entry(seqid, source, feature, int(start), int(end), score, strand, phase, attr)
#edef

###############################################################################

class GFF3(object):
  entries = None
  fileName = None
  seqids = None
  index = None
  topLevel = None
  idField = None
  parentField = None
  nameField = None

  def __init__(self, entries=None, fileName=None, idField="ID", parentField="Parent", nameField="ID", fileSkipLines=0, fileMaxLines=None):
    if (fileName is not None) and entries is None:
      self.fileName = fileName
      entries = readGFF3File(fileName, skipLines=fileSkipLines, maxLines = fileMaxLines)
    #fi

    self.idField = idField
    self.parentField = parentField
    self.nameField = nameField

    if entries is not None:
      self.entries = entries
      self.seqids  = set([ e.seqid for e in self.entries])
      self.index, self.topLevel = self.__index()
    #fi
  #edef

  def __str__(self):
    dstr  = "GFF3 object\n"
    dstr += " Where: %s\n" % self.fileName
    dstr += " Entries: %d\n" % len(self.entries)
    dstr += " Top level statistics:\n"
    for featureType in self.topLevel:
      dstr += "  * %s : %d\n" % (featureType, len(self.topLevel[featureType]))
    #efor
    return dstr
  #edef
      

  def __index(self):
    idx = {}
    topLevel = {}
    for i, e in enumerate(self.entries):
      ID = e.attr[self.idField] if self.idField in e.attr else None
      if ID is None:
        continue
      #fi

      parent = e.attr[self.parentField] if self.parentField in e.attr else None
      if ID not in idx:
        idx[ID] = [i, [] ]
      else:
        idx[ID] = [i, idx[ID][1] ]
      #fi
      if parent is None:
        if e.feature not in topLevel:
          topLevel[e.feature] = []
        #fi
        topLevel[e.feature].append(ID)
      else:
        if parent not in idx:
          idx[parent] = [ None, [] ]
        #fi
        idx[parent][1].append((i, ID))
      #fi
    #efor
    return idx, topLevel
  #edef

  def getIDEntry(self, ID):
    if ID in self.index:
      return self.entries[self.index[ID][0]]
    else:
      return None
    #fi
  #edef

  def getIDChildren(self, ID):
    if ID in self.index:
      return self.index[ID][1]
    else:
      return None
    #fi
  #edef
      
  def getChildren(self, ID, feature=None, depth=None, containParent=False, newObject=False):

    relEntries = [ self.getIDEntry(ID) ] if containParent else []
    ids = [ (0, c[0], c[1]) for c in self.getIDChildren(ID) ]

    while len(ids) > 0:
      cdepth, cidIndex, cid = ids.pop(0)
      relEntries.append(self.entries[cidIndex])
      if (cid in self.index) and ( (depth is None) or (cdepth+1 < depth)):
        ids.extend([ (cdepth+1, c[0], c[1]) for c in self.getIDChildren(cid) ])
      #fi
    #ewhile

    if newObject:
      return GFF3(entries = relEntries)
    else:
      return relEntries
    #fi
  #edef  

  def indexByInterval(self):
    try: 
      from intervaltree import Interval, IntervalTree
    
      t = { seqid: IntervalTree() for seqid in self.seqids }
      for e in [ e for e in self.entries if e.type.lower() == "mrna"]:
        t[e.seqid][e.start:e.end] = e
      #efor
      return t
    except ImportError:
      return {}
  #edef

  def getID(self, ID):
    if ID in self.index:
      return self.entries[self.index[ID][0]]
    else:
      return None
    #fi
  #edef

  def areTandem(self, id1, id2):
    f1 = self.getID(id1)
    f2 = self.getID(id2)
    if gene1.seqid != gene2.seqid:
      return False
    #fi

    inregion = set([ e[-1].attr["ID"] for e in self.interval[gene1.seqid][min(gene1.end,gene2.end):max(gene1.start,gene2.start)] ])
    if len(inregion - set([gene1.attr["ID"], gene2.attr["ID"]])) == 0:
      return True
    else:
      return False
    #fi
  #edef

  def areSameStrand(self, id1, id2):
    f1 = self.getID(id1)
    f2 = self.getID(id2)
    return f1.strand == f2.strand
  #edef

  def write(self, fileName):
    with (gzip.open(fileName, "wb") if fileName[-2:] == "gz" else open(fileName, "w")) as gffFile:
      gffFile.write("#gff-version\t3\n")
      for e in self.entries:
        #gff3Entry = namedtuple("gff3Entry", "seqid, source, feature, start, end, score, strand, phase, attr")
        gffFile.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (e.seqid, e.source, e.feature, e.start, e.end, e.score, e.strand, e.phase, ';'.join([ '%s=%s' % (k,v) if (v is not None) else k for (k,v) in e.attr.items() ]) ))

#eclass

###############################################################################

def readGFF3File(filename, skipLines=0, maxLines=None):
  G = []
  nLines = 0
  with (gzip.open(filename, "rt") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      nLines += 1
      if nLines < skipLines:
        continue
      elif (maxLines is not None) and nLines > maxLines:
        break
      #fi

      if len(row) != 9:
        continue
      #fi
      G.append(parseGFF3entry(row))
    #efor
  #ewith
  return G
#edef


###############################################################################
