from .. import GFF3
from .. import Fasta
from ..lazyObject import LazyObject
from .. import utils as utils
import pandas as pd
import numpy as np
import tabix
import vcf
from collections import namedtuple

###############################################################################

class ResourceManager(LazyObject):

  _resource = None
  _initialized = None
  _fmObject = None
  _requiredFiles = None

  def __init__(self, fmObject, requiredFiles):
    print("Initializing the class only NOW")
    self._fmObject = fmObject
    self._requiredFiles = requiredFiles
    if not(fmObject.satisfyRequiredFiles(requiredFiles)):
      print("Error: could not initialize")
      self._initialized = False
    else:
      self._initialized = True
    #fi
  #edef

  def __call__(self):
    return self._resource
  #edef

  def __str__(self):
    dstr =  "ResourceManager Object:\n"
    dstr += " Initialized: %s\n" % ("Yes" if self._initialized else "No")
    dstr += " Resources managed: \n"
    for fileID in self._requiredFiles:
      dstr += "  * %s -> %s\n" % (fileID, self._fmObject.getFileName(fileID))
    #efor
    return dstr
  #edef

#eclass

###############################################################################

class TSVResourceManager(ResourceManager):
  def __init__(self, fmObject, tsvFile, fieldNames=None, delimiter='\t', skiprows=1):
    ResourceManager.__init__(self, fmObject, [ tsvFile ])
    if self._initialized:
      self._resource = pd.read_csv(self._fmObject.getFileName(tsvFile), delimiter='\t', names=fieldNames, skiprows=1)
    #fi
  #edef
#eclass

###############################################################################

class TabixTSVResourceManager(ResourceManager):

  def __init__(self, fmObject, tsvFile, tabixFile, fieldNames=None):
    ResourceManager.__init__(self, fmObject, [ tsvFile, tabixFile ])
    if self._initialized:
      self._resource = tabix.open(self._fmObject.getFileName(tsvFile))
      if fieldNames is not None:
        self._namedtuple = namedtuple("tabixTsv", fieldNames)
      #fi
      self._fieldNames = fieldNames
    #fi
  #edef

  def query(self, seqid, start, end, pandas=False, namedtuple=False, tabixIter=True):
    res = utils.tabixQueryWrapper(self._resource, seqid, start, end)
    if pandas:
      return pd.DataFrame([ r for r in res])
    elif namedtuple and (self._fieldNames is not None):
      return [ self._namedtuple(*r) for r in res ]
    elif tabixIter:
      return res
    else:
      return [ r for r in res ]
    #fi
  #edef
#eclass

###############################################################################

class VCFResourceManager(ResourceManager):
  def __init__(self, fmObject, vcfFile, tabixFile, fieldNames=None):
    ResourceManager.__init__(self, fmObject, [ vcfFile, tabixFile ])
    if self._initialized:
      self._resource = vcf.Reader(filename=self._fmObject.getFileName(vcfFile), compressed=True)
    #fi
  #edef

  def query(self, seqid, start, end, pandas=False, namedtuple=False, tabixIter=True):
    res = utils.vcfQueryWrapper(self._resource, seqid, start, end)
    return res
  #edef
#eclass

###############################################################################

class GFF3ResourceManager(ResourceManager, GFF3):
  def __init__(self, fmObject, gff3File, **kwargs):
    ResourceManager.__init__(self, fmObject, [ gff3File ])
    if self._initialized:
      GFF3.__init__(self, fileName=self._fmObject.getFileName(gff3File), **kwargs)
    #fi
  #edef

  def __str__(self):
    return GFF3.__str__(self)
  #edef
#eclass

###############################################################################

class FastaResourceManager(ResourceManager, Fasta):
  def __init__(self, fmObject, fastaFile, fileSkipLines=0, fileMaxLines=None):
    ResourceManager.__init__(self, fmObject, [ fastaFile ])
    if self._initialized:
      Fasta.__init__(self, fileName=self._fmObject.getFileName(fastaFile))
    #fi
  #edef

  def __str__(self):
    return Fasta.__str__(self)
  #edef
#eclass

###############################################################################
