from ..structures import Dataset
from .. import formats
from .. import utils
from .. import stats
from ..config import settings

pd = utils.py.loadExternalModule("pandas")

import os

###############################################################################

# https://www.biostars.org/p/71737/


class KEGG(Dataset):

  versions = {
    "human" : "hsa",
    "mouse" : "mmu",
    "drosophilia" : "dme"
  }

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fileIndex = self.__genFileIndex(version)
    Dataset.__init__(self, fileIndex)
    self.version = version
    self.__orgID = self.versions[self.version]

    self._registerObject("_orgMap", formats.TSVMap, [ "org_map"], fileIndex["org_map"].path, delimiter='\t')
    self._registerObject("_featureData", formats.SQLDict, ["feature_data"], fileIndex["feature_data"].path)

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef


  def __genFileIndex(self, version, where=None):
     finalPath = '%s/kegg/%s' % ( (settings.getDataDir() if where is None else where), version)
     orgKey = self.versions[version]
     files = {}
     files['org_map'] = utils.Acquire(where=where).curl("http://rest.kegg.jp/link/%s/pathway" % orgKey).finalize('%s/org_map.tsv' % finalPath)
     files['feature_data'] = utils.Acquire(where=where).touch('%s/feature_data.dict.sqlite' % finalPath)
     return files
  #edef

  ###############################################################################

  def getPathways(self):
    return list(self._orgMap.fromKeys)
  #edef

  def getGenes(self):
    return list(self._orgMap.toKeys)
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
        return "path:%s%05d" % (self.__orgID, ID)
      elif ID.isdigit():
        return "path:%s%05d" % (self.__orgID, int(ID))
      elif ID[:5+len(self.__orgID)] == "path:%s" % self.__orgID:
        return ID
      elif ID[:len(self.__orgID)] == self.__orgID:
        return "path:%s" % ID
      else:
        utils.dbm("Don't know what to do with: %s" % str(ID))
        return ""
      #fi
    else: # We want a gene ID
      if isinstance(ID, int):
        return "%s:%d" % (self.__orgID, ID)
      elif ID.isdigit():
        return "%s:%d" % (self.__orgID, int(ID)) # Remove leading zeros from string
      elif ID[:1+len(self.__orgID)] == '%s:' % self.__orgID:
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
      url = "http://rest.kegg.jp/get/%s" % (ID)
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
