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
    "drosophilia" : "dme",
    "yeast" : 'sce'
  }

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    """
      Generate a KEGG dataset

      Inputs:
       - version : Which organism to use
       - **kwargs : Dataset structure arguments
      Output:
       KEGG data structure
    """

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
     orgKey = self.versions[version] if version in self.versions else version
     files = {}
     files['org_map'] = utils.Acquire(where=where).curl("http://rest.kegg.jp/link/%s/pathway" % orgKey).finalize('%s/org_map.tsv' % finalPath)
     files['feature_data'] = utils.Acquire(where=where).touch('%s/feature_data.dict.sqlite' % finalPath)
     return files
  #edef

  ###############################################################################

  def getPathways(self):
    """ Get a list of all pathway IDs """
    return list(self._orgMap.fromKeys)
  #edef

  def getGenes(self):
    """ Get all Kegg gene IDs """
    return list(self._orgMap.toKeys)
  #edef

  def getGeneIDs(self):
    """ Return the entrez GeneIDs from the kegg map, rather than the kegg IDs (with hsa: infront)"""
    return [ g.split(':')[1] for g in self.getGenes() ]
  #edef

  def getPathwayGenes(self, pathwayID):
    """
      getPathwayGenes: Get all kegg gene IDs for a given pathway

      Inputs:
       - pathwayID : Kegg Pathway ID
      Outputs:
        List of kegg gene IDs
    """
    return self._orgMap.lookup(self._formatFeatureID(pathwayID, True))
  #edef

  def getPathwayGeneIDs(self, pathwayID):
    """
      getPathwayGeneIDs: Get all entrez gene IDs for a given pathway

      Inputs:
       - pathwayID : Kegg Pathway ID
      Outputs:
        List of entrez gene IDs
    """
    return [ g.split(':')[1] for g in self.getPathwayGenes(pathwayID) ]
  #edef

  def getGenePathways(self, geneID):
    """
      getGenePathways: Get all pathway IDs for a specific gene

      Inputs:
       - geneID : entrez GeneID, or KEGG gene ID
      Outputs:
        List of pathway IDs
    """
    return self._orgMap.inverse(self._formatFeatureID(geneID, False))
  #edef

  def _formatFeatureID(self, ID, pathway):
    if ID is None:
      raise ValueError("None is not a valid ID for KEGG IDs.")
    #fi
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
      elif ID in self.getGeneIDs():
        return "%s:%s" % (self.__orgID, ID)
      else:
        utils.dbm("Don't know what to do with '%s'" % str(ID))
      #fi
    #fi
  #edef

  def getPathwayInfo(self, pathwayID):
    """
      getPathwayInfo: get Information of a pathway

      Inputs:
        pathwayID: KEGG pathway identifier
      Outputs:
        Text pathway information
    """
    return self._getFeature(self._formatFeatureID(pathwayID, True))
  #edef

  def getPathwayName(self, pathwayID):
    """
      getPathwayName: get name of a pathway

      Inputs:
        pathwayID: KEGG pathway identifier
      Outputs:
        Text pathway name
    """
    return self.getPathwayInfo(pathwayID).split('\n')[1][4:].strip()
  #edef

  def getGeneInfo(self, geneID):
    """
      getGeneInfo: get Information of a gene

      Inputs:
        geneID : entrez geneID, or kegg gene identifier
      Outputs:
        Text gene information
    """
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

  def enrich(self, yourSet, background=None, pathway=None, abcd_values=False, method=None, **kwargs):
    """
    Enrich: Check enrichment of KEGG pathways in a given set

    Inputs:
      - yourSet: List of Entrez Gene IDs to test
      - pathway: List of pathways (or single pathway) to test (Defaults to all pathways that your geneIDs are present in)
      - background: List of Entrez Gene IDs to use as background (e.g. set of all expressed genes)
      - abcd_values: Boolean. If True, it will return the actual element values in the contingency table, rather than just counts
      - method: Type of multiple testing correction procedure to use
      - **kwargs: Additional arguments for biu.stats.p_adjust

    Outputs:
     - df : Pandas Data Frame of test results
    """
    yourSet = set([ str(ID) for ID in yourSet if ID is not None]) & set(self.getGeneIDs())

    if pathway is None:
        pathway = set([])
        for gene in yourSet:
          pathway = set(self.getGenePathways(gene)) | pathway
        #efor
    #fi
    
    if background is None:
        background = self.getGeneIDs()
    #fi
    background = set(background)
    
    if isinstance(pathway, str):
        pathway = [ pathway ]
    #fi
    R = []
    B = self.getGeneIDs()
    for p in pathway:
        pathwayGenes = set(self.getPathwayGeneIDs(p)) & background
        res = stats.enrichment.setEnrichment(yourSet, pathwayGenes, B, abcd_values=abcd_values)
        R.append((p, res.method, res.c2statistic, res.oddsratio, res.pvalue, res.table))
    #efor

    df = pd.DataFrame(R, columns=['pathway', 'method', 'c2statistic', 'oddsratio', 'p', 'table'])
    if method is not None:
      df['q'] = stats.p_adjust(df.p.values, method, **kwargs)
    #fi

    return df
  #edef

#eclass
