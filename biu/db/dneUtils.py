from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

import os

###############################################################################

def urlFileIndex():
  files = {}

  files["probs"]       = ("https://media.nature.com/original/nature-assets/ng/journal/v46/n9/extref/ng.3050-S2.xls", "probs.xls", {})
  files["constrained"] = ("https://media.nature.com/original/nature-assets/ng/journal/v46/n9/extref/ng.3050-S3.xls", "constrained.xls", {})

  return { k : (u, 'DNE/%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class DNE(fm.FileManager):

  _probs = None
  _constrained = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=[ "_probs", "_constrained" ], skiprows=0, **kwargs)

    self._probs = rm.XLSResourceManager(self, "probs")
    self._constrained = rm.XLSResourceManager(self, "constrained")
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

    probs = self._probs["mutation_probabilities"]
    probs = probs[probs.gene.apply(filtFunc)]

    constrained = self._constrained["constrained_genes"]
    constrained = constrained[constrained.gene.apply(filtFunc)]

    return probs.join(constrained.set_index("transcript").drop(["gene", "bp"], axis=1), how="left", on="transcript")
  #edef

  def __contains__(self, key):
    return (key in set(self._constrained["constrained_genes"].gene))
  #edef

#eclass
