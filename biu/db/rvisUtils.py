from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils

import os

###############################################################################

# https://www.biostars.org/p/71737/

def urlFileIndex():
  files = {}

  files["original"] = ("http://genic-intolerance.org/data/GenicIntolerance_v3_12Mar16.txt", "rvis_v9.tsv", {})
  files["exAc"]     = ("http://genic-intolerance.org/data/RVIS_Unpublished_ExAC_May2015.txt", "rvis_exac.tsv", {})
  files["exAc2"]    = ("http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt", "rvis_exac2.tsv", {})

  return { k : (u, 'rvis/%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class RVIS(fm.FileManager):

  _original = None
  _exAc     = None
  _exAc2    = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=[ "_original", "_exAc", "_exAc2" ], skiprows=0, **kwargs)

    self._original = rm.TSVResourceManager(self, "original", delimiter='\t')
    self._exAc     = rm.TSVResourceManager(self, "exAc", delimiter='\t')
    self._exAc2    = rm.TSVResourceManager(self, "exAc2", delimiter='\t')
  #edef

  __getitem__colMap = {
    "CCDS_r9"   : "GENE",
    "ALL_0.1%"  : "RVIS",
    "%ALL_0.1%" : "RVIS_p",
    "%ExAC_0.05%popn" : "RVIS_exAc_p",
    "LoF-FDR[ExAC]"   : "LoF_p",
    "CCDSr20" : "GENE",
    "RVIS[pop_maf_0.05%(any)]" : "RVIS_gnomad",
    "%RVIS[pop_maf_0.05%(any)]" : "RVIS_gnomad_p",
    "Edge_case_RVIS[pop_maf_0.05%(any)]" : "edge_case",
    "%OE-ratio_[ExAC v2]"                : "edge_case_oe_ratio"
  }

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

    # RVIS     %ExAC RVIS     ExAC LoF FDR     %ExAC v2 RVIS     Edge Case (%OE-ratio) 

    relO = self._original[self._original.GENE.apply(filtFunc)][["GENE", "ALL_0.1%", "%ALL_0.1%", "%ExAC_0.05%popn", "LoF-FDR[ExAC]"]]
    relG = self._exAc2[self._exAc2.CCDSr20.apply(filtFunc)][["CCDSr20", "RVIS[pop_maf_0.05%(any)]", "%RVIS[pop_maf_0.05%(any)]", "Edge_case_RVIS[pop_maf_0.05%(any)]", "%OE-ratio_[ExAC v2]" ]]
    
    relO = relO.rename(columns=self.__getitem__colMap)
    relG = relG.rename(columns=self.__getitem__colMap)
    
    return relO.join(relG.set_index("GENE"), how="left", on="GENE")
  #edef

#eclass

