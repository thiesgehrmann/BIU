from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from .. import utils
from .. import stats

import pandas as pd

import os

###############################################################################

# https://www.biostars.org/p/71737/

versions = {
  "human" : "hsa",
  "mouse" : "mmu",
  "drosophilia" : "dme"
}

def urlFileIndex(version):
  files = {}

  files["org_map"] = ( "http://rest.kegg.jp/link/%s/pathway" % versions[version], "org_map.tsv", {})
  files["feature_data"] = ("http://rest.kegg.jp/get/%s", "feature_data.sqlite", {})

  return { k : (u, 'kegg_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for version in versions:
    print(" * %s" % version)
  #efor
#edef

###############################################################################

class KEGG(fm.FileManager):

  _orgMap = None
  _featureData = None
  _orgID = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=[ "_orgMap", "_featureData" ], skiprows=0, **kwargs)
    self.version = version
    self._orgID = versions[self.version]

    self._orgMap = rm.TSVMapResourceManager(self, "org_map", delimiter='\t')
    self._featureData = rm.SQLDictResourceManager(self, "feature_data")

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  ###############################################################################

  def getPathways(self):
    return list(self._orgMap.lookupKeys())
  #edef

  def getGenes(self):
    return list(self._orgMap.inverseKeys())
  #edef

  def getGeneIDs(self):
    """ Return the NCBI GeneIDs from the kegg map, rather than the kegg IDs (with hsa: infront)"""
    return [ g.split(':')[1] for g in self.getGenes() ]
  #edef

  def getPathwayGenes(self, pathwayID):
    return self._orgMap.lookup(self._formatFeatureID(pathwayID, True))
  #edef

  def getPathwayGeneIDs(self, pathwayID):
    return [ g.split(':')[1] for g in self.getPathwayGenes(pathwayID) ]
  #edef

  def getGenePathways(self, geneID):
    return self._orgMap.inverse(self._formatFeatureID(geneID, False))
  #edef

  def _formatFeatureID(self, ID, pathway):
    if pathway:
      if isinstance(ID, int):
        return "path:%s%05d" % (self._orgID, ID)
      elif ID.isdigit():
        return "path:%s%05d" % (self._orgID, int(ID))
      elif ID[:5+len(self._orgID)] == "path:%s" % self._orgID:
        return ID
      elif ID[:len(self._orgID)] == self._orgID:
        return "path:%s" % ID
      else:
        utils.dbm("Don't know what to do with: %s" % str(ID))
        return ""
      #fi
    else: # We want a gene ID
      if isinstance(ID, int):
        return "%s:%d" % (self._orgID, ID)
      elif ID.isdigit():
        return "%s:%d" % (self._orgID, int(ID)) # Remove leading zeros from string
      elif ID[:1+len(self._orgID)] == '%s:' % self._orgID:
        return ID
      else:
        utils.dbm("Don't know what to do with '%s'" % str(ID))
      #fi
    #fi
  #edef

  def getPathwayInfo(self, pathwayID):
    return self._getFeature(self._formatFeatureID(pathwayID, True))
  #edef

  def getPathwayName(self, pathwayID):
    return self.getPathwayInfo(pathwayID).split('\n')[1][4:].strip()
  #edef

  def getGeneInfo(self, geneID):
    return self._getFeature(self._formatFeatureID(geneID, False))
  #edef

  def _getFeature(self, ID):
    if ID in self._featureData:
      return self._featureData[ID]
    else:
      url = self._fileIndex["feature_data"][0] % (ID)
      utils.dbm("Downloading via REST from '%s'" % url)
      dat = str(utils.getCommandOutput("curl --silent -L '%s'" % url).decode('UTF-8'))
      self._featureData[ID] = dat
      return dat
    #fi
  #edef

  def enrich(self, yourSet, pathway=None, correctionType=None, **kwargs):
    if pathway is None:
        pathway = set([])
        for gene in yourSet:
          pathway = set(self.getGenePathways(gene)) | pathway
        #efor
    #fi
    if isinstance(pathway, str):
        pathway = [ pathway ]
    #fi
    R = []
    B = self.getGeneIDs()
    for p in pathway:
        pathwayGenes = self.getPathwayGeneIDs(p)
        res = stats.enrichment.setEnrichment(yourSet, pathwayGenes, B)
        R.append((p, res.method, res.c2statistic, res.oddsratio, res.pvalue))
    #efor

    df = pd.DataFrame(R, columns=['pathway', 'method', 'c2statistic', 'oddsratio', 'p'])
    if correctionType is not None:
      df['q'] = stats.correction.correct(df.p.values, correctionType, **kwargs)
    #fi

    return df
#edef

#eclass
