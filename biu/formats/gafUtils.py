from .. import utils

import pandas as pd

class GAF(object):

  # GO annotation File: http://geneontology.org/page/go-annotation-file-gaf-format-21

  _fieldNames = [ "db", "db_o_id", "db_o_symbol", "qualifier", "go_id", "db_ref", "evidence", "w_o_f", "aspect", "db_o_name", "db_o_synonym", "db_o_type", "date", "assigned_by", "annot_ext", "gene_product_form_id" ]

  _entries = None
  _fileName = None
  _annotIndex = None
  _objectIndex = None

  def __init__(self, fileName, **kwargs):
    self._entries = pd.read_csv(fileName, index_col=False, names=self._fieldNames, **kwargs)
    self._fileName = fileName
    self._annotIndex = {}
    self._objectIndex = {}

    for row in self._entries.itertuples():
      i = int(row.Index)
      goID = row.go_id
      objectID = row.db_o_id

      if goID not in self._annotIndex:
        self._annotIndex[goID] = []
      #fi
      if objectID not in self._objectIndex:
        self._objectIndex[objectID] = []
      #fi

      self._annotIndex[goID].append( (objectID, i) )
      self._objectIndex[objectID].append( (goID, i) )
    #efor
  #edef

  def getAnnots(self, objectID):
    if objectID not in self._objectIndex:
      return []
    else:
      return [ a[0] for a in self._objectIndex[objectID] ]
  #edef

  def getAnnotated(self, annotID):
    if annotID not in self._annotIndex:
      return []
    else:
      return [ a[0] for a in self._annotIndex[annotID] ]
    #fi
  #edef

  def __str__(self):
    dstr  = "GAF (GO Annotation File) Object\n"
    dstr += " Where: %s\n" % self._fileName
    dstr += " # Annotations : %d\n" % self._entries.shape[0]
    dstr += " # Objects     : %d\n" % len(self.objectIndex.keys())
    dstr += " # GO terms    : %d\n" % len(self.annotIndex.keys())
    return dstr
  #edef

#eclass
