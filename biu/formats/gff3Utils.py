from .. import utils

import errno
import os
import csv
import gzip
from collections import namedtuple

import intervaltree

import gzip

###############################################################################

idField     = lambda e: e.attr["ID"] if "ID" in e.attr else e.attr["Name"] if "Name" in e.attr else None
parentField = lambda e: e.attr["Parent"] if "Parent" in e.attr else None
nameField   = lambda e: e.attr["Name"] if "Name" in e.attr else None
seqIDField  = lambda e: e.seqid

###############################################################################

class GFF3Entry(object):

  __slots__ = [ 'seqid', 'source', 'feature', 'start', 'end', 'score', 'phase', 'strand', 'attr', '__idField', '__parentField', '__nameField' ]

  #seqid = None
  #source = None
  #feature = None
  #start = None
  #end = None
  #score = None
  #phase = None
  #attr = None

  def __init__(self, row, idField=idField, parentField=parentField, nameField=nameField, **kwargs):
  
    self.__idField = idField
    self.__parentField = parentField
    self.__nameField = nameField

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
  
    (self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, attr) = row
    if isinstance(attr, dict):
      self.attr = attr
    else:
      self.attr = dict([attrsplit(x.strip()) for x in attr.split(";") ])
    #fi
    self.start = int(self.start)
    self.end   = int(self.end)
  #edef

  @property
  def id(self):
    return self.__idField(self)
  #edef

  @property
  def parent(self):
    return self.__parentField(self)
  #edef

  @property
  def name(self):
    return self.__nameField(self)
  #edef

  def seq(self, fastaObject):
    if self.seqid not in fastaObject:
      utils.error("Sequence '%s' not found in provided fastaObject." % self.seqid)
      return None
    #fi

    substr = fastaObject[self.seqid][self.start-1:self.end]
    if self.strand == '-':
      substr = substr.revcomp()
    #fi

    return substr
  #edef

  def __attrString(self):
    return ';'.join([ '%s=%s' % (k,v) if (v is not None) else k for (k,v) in self.attr.items() ])

  def outputString(self):
    return "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s" % (self.seqid, self.source, self.feature, self.start, 
                                                   self.end, self.score, self.strand, self.phase,
                                                   self.__attrString() )
  #edef

  def copy(self):
    return GFF3Entry([self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, self.attr], idField=self.__idField, parentField=self.__parentField, nameField=self.__nameField)
  #edef

  def withoutParent(self):
    newAttr = { k : v for (k,v) in self.attr.items() if v != self.parent }
    return GFF3Entry([self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, newAttr], idField=self.__idField, parentField=self.__parentField, nameField=self.__nameField)
  #edef

  def __str__(self):
    dstr = "GFF3Entry(seqid:%s, source:%s, feature:%s, start:%d, end:%d, score:%s, strand:%s, phase:%s, attr:%s)" % (self.seqid, self.source, self.feature, self.start,
                                                                                                                     self.end, self.score, self.strand, self.phase,
                                                                                                                     self.__attrString() )
    return dstr
  #edef

  def __repr__(self):
    return self.__str__()
#eclass

###############################################################################

class GFF3(object):

  __slots__ = [ 'entries', 'seqids', 'index', 'features', '__index', '__intervalIndex', '__fileName' ]

  #entries = None
  #seqids = None
  #index = None
  #features = None

  #__index = None
  #__intervalIndex = None
  #__fileName = None

  def __init__(self, data, **kwargs):

    if isinstance(data, str):
      utils.dbm("GFF input source is file.")
      self.__fileName = data
      self.entries = GFF3.read(data, **kwargs)
    elif isinstance(data, type(self)):
      utils.dbm("GFF input source is GFF3 structure")
      self.entries = data.entries
    else:
      utils.dbm("GFF input source is list of GFF3Entries.")
      self.entries = data
    #fi
    self.seqids  = set([ e.seqid for e in self.entries])
    self.__index, self.features = self._index()
  #edef

  def __iter__(self):
    return self.entries.__iter__()
  #edef

  def __str__(self):
    dstr  = "GFF3 object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    dstr += " Entries: %d\n" % len(self.entries)
    dstr += " Indexed: %s\n" % ("Yes" if self.__index is not None else "No")

    if self.__intervalIndex is not None:
      dstr += " Interval Indexes:\n"
      for indexFeatures in self.__intervalIndex:
        dstr += "  * %s\n" % ",".join(indexFeatures)
      #efor
    #fi

    dstr += " Feature statistics:\n"
    for featureType in self.features:
      dstr += "  * %s : %d\n" % (featureType, len(self.features[featureType]))
    #efor
    return dstr
  #edef

  def __getitem__(self, c):
    if isinstance(c, str):
      return self.getIDEntry(c)
    else:
      return self.entries[c]
    #fi
  #edef

  def _index(self):
    internal_counter = 0

    idx = {}
    features = {}
    for i, e in enumerate(self.entries):
      ID = e.id
      if ID is None:
        ID = "internal.%d" % internal_counter
        internal_counter += 1
      #fi

      if ID not in idx:
        idx[ID] = [i, [] ]
      else:
        idx[ID] = [i, idx[ID][1] ]
      #fi

      # Add feature to top Level index
      feature = e.feature
      if feature != "":
        if feature not in features:
          features[feature] = []
        #fi
        features[e.feature].append(ID)
      #fi

      # Construct hierarchical structure
      parent = e.parent
      if parent is not None:
        if parent not in idx:
          idx[parent] = [ None, [] ]
        #fi
        idx[parent][1].append((i, ID))
      #fi
    #efor
    return idx, features
  #edef

  def __len__(self):
    return len(self.entries)
  #edef

  def seq(self, ID, fastaObject, feature='CDS'):
    entries = self.getChildren(ID, feature=feature).entries
    if len(entries) == 0:
      utils.error("Could not find ID '%s'" % ID)
      return None
    #fi
    entries = sorted(entries, key=lambda e: e.start)
    if entries[0].strand == '-':
      entries = entries[::-1]
    #fi

    sequence = sum([ e.seq(fastaObject) for e in entries ])
    sequence.setName(ID)
    return sequence
  #edef

  def getIDEntry(self, ID):
    if ID in self.__index:
      return self.entries[self.__index[ID][0]]
    else:
      return None
    #fi
  #edef

  def getIDChildren(self, ID):
    if ID in self.__index:
      return self.__index[ID][1]
    else:
      utils.error("No children found for '%s'" % ID)
      return []
    #fi
  #edef
      
  def getChildren(self, ID, feature=None, depth=None, containParent=False):

    relEntries = [ self.getIDEntry(ID).withoutParent() ] if containParent else []
    ids = [ (0, c[0], c[1]) for c in self.getIDChildren(ID) ]

    while len(ids) > 0:
      cdepth, cidIndex, cid = ids.pop(0)
      relEntries.append(self.entries[cidIndex])
      if (cid in self.__index) and ( (depth is None) or (cdepth+1 < depth)):
        ids.extend([ (cdepth+1, c[0], c[1]) for c in self.getIDChildren(cid) ])
      #fi
    #ewhile

    if feature is not None:
      relEntries = [ r for r in relEntries if r.feature in feature ]
    #fi

    # Strip the parent tag if the parent == ID
    # "seqid, source, feature, start, end, score, strand, phase, attr"
    if not(containParent):
      E = []
      for e in relEntries:
        if e.parent == ID:
          e = e.withoutParent()
          E.append(e)
        else:
          E.append(e)
        #fi
      #efor
      relEntries = E
    #fi

    return GFF3(data = relEntries)
  #edef  

  def __getIntervalIndex(self, features):
    if features is None:
      features = list(self.features.keys())
    #fi

    if isinstance(features, str):
      features = [ features ]
    #fi
    features = tuple(sorted([ f.lower() for f in features]))

    if self.__intervalIndex is None:
      self.__intervalIndex = {}
    #fi

    if features not in self.__intervalIndex:
      try: 
        from intervaltree import Interval, IntervalTree
      
        t = { str(seqid) : IntervalTree() for seqid in self.seqids }
        for e in [ e for e in self.entries if e.feature.lower() in features ]:
          t[str(e.seqid)][e.start:e.end] = e
        #efor
      except ImportError:
        return {}
      #etry
      self.__intervalIndex[features] = t
    #fi
    return self.__intervalIndex[features]
  #edef

  def query(self, chromosome, start, stop, features=None):
    return self.queryRegions([(chromosome, start, stop)], features=features)
  #edef

  def queryRegions(self, regions, features=None):
    R = []
    for (c,s,e) in regions:
      r = self.__query(c, s, e, features)
      R.extend(r)
    #efor
    return GFF3(R)
  #edef

  def __query(self, c, s, e, features):
    # Make a string of the seqid
    c = str(c)
    # reorder the coords
    if s > e:
      s, e = ( min(s, e), max(s, e))
    #fi

    relIntervalIndex = self.__getIntervalIndex(features)
    if c not in relIntervalIndex:
      utils.warning("Seqid '%s' not in GFF." % c)
      return []
    #fi
    return [ i[-1] for i in relIntervalIndex[c][s:e] ]
  #edef

  def getID(self, ID):
    if ID in self.__index:
      return self.entries[self.__index[ID][0]]
    else:
      return None
    #fi
  #edef

  @property
  def dataFrame(self):
    import pandas as pd
    return pd.DataFrame([ (e.seqid, e.source, e.feature, e.start, e.end, e.score, e.phase, e.attr) for e in self.entries ],
                        columns=("seqid", "source", "feature", "start", "end", "score", "phase", "attr") )
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
        gffFile.write("%s\n" % e.outputString())
      #efor
    #ewith
  #edef

###############################################################################

  @staticmethod
  def read(filename, skipLines=0, maxLines=None, allowAdditionalColumns=False, **kwargs):
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
  
        # Unless we allow it with allowAdditionalColumns, we require exactly 9 columns.
        ncolumns = len(row)
        if (ncolumns < 9) or ((ncolumns > 9) and not(allowAdditionalColumns)):
          continue
        #fi
        G.append(GFF3Entry(row[:9], **kwargs))
      #efor
    #ewith
    return G
  #edef

#eclass
###############################################################################
