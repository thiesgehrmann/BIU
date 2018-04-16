from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

import os

###############################################################################

def urlFileIndex():
  files = {}

  files["gdi"] = ("http://lab.rockefeller.edu/casanova/assets/file/GDI_full_10282015.txt", "gdi.tsv", {})

  return { k : (u, 'GDI/%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class GDI(fm.FileManager):

  _gdi = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=[ "_gdi" ], skiprows=0, **kwargs)

    self._gdi = rm.TSVResourceManager(self, "gdi", delimiter='\t')
  #edef

  def __getitem__(self, geneSymbols):

    if geneSymbols is None:
      filtFunc = lambda g: True
    else:
      if isinstance(geneSymbols, str):
        geneSymbols = [ geneSymbols ]
      #fi
      geneSymbols = set(geneSymbols)
      filtFunc = lambda g: g in geneSymbols
    #fi

    return self._gdi[self._gdi.Gene.apply(filtFunc)]
  #edef

#eclass
