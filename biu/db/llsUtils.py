from ..structures import Dataset
from ..config import settings as settings
from .. import formats
from .. import utils

import itertools

pd = utils.py.loadExternalModule("pandas")

###############################################################################

class LLS(Dataset):

  versions = { "current":
    { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "M", "X", "Y" ] }
  }
  version = None
  vcf = None

  def __init__(self, version=list(versions.keys())[0], where="/exports/molepi/LLSSEQ", localCopy={}):
    fileIndex = self.__genFileIndex(version, where)

    Dataset.__init__(self, fileIndex, localCopy=localCopy)
    self.version = version

    for chrID in self.versions[self.version]["chrs"]:
      self._registerObject('vcf_%s' % chrID, formats.VCF, [ "vcf_%s" % chrID, "vcf_%s_tbi" % chrID ], fileIndex["vcf_%s" % chrID].path, tabix=True)
    #efor

    self._registerObject('percentiles', pd.read_csv, [ "percentiles" ], fileIndex["percentiles"].path, sep='\t')
    self._registerObject('ids', pd.read_csv, ["ids"], fileIndex["ids"].path)
    self._registerObject('phen', pd.read_csv, [ "phen"], fileIndex["phen"].path, delimiter=' ')

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     files = {}
     for chrID in self.versions[version]["chrs"]:
       files['vcf_%s' % chrID] = utils.Acquire("%s/tbx/merged.chr%s.vcf.bgz" % (where, chrID), where=where)
       files['vcf_%s_tbi' % chrID] = utils.Acquire("%s/tbx/merged.chr%s.vcf.bgz.tbi" % (where, chrID), where=where)
     #efor

     files['phen'] = utils.Acquire("/exports/molepi/tgehrmann/data/LLS_data/lls_phenotypes218.txt")
     files['percentiles'] = utils.Acquire("/exports/molepi/tgehrmann/data/LLS_data/lls_percentiles.tsv")
     files['ids'] = utils.Acquire('/exports/molepi/tgehrmann/data/LLS_data/lls_ids.csv')

     return files
  #edef

  def query(self, chrID, start, end, **kwargs):
    chrID = str(chrID)
    oname = "vcf_%s" % chrID
    if self._objectExists(oname):
      return self._getObject(oname).query(chrID, start, end, **kwargs)
    else:
      utils.error("Could not find chromosome '%s'" % chrID)
      return iter(())
    #fi
  #edef

  def queryRegions(self, regions, extract=None, **kwargs):
    R = []
    template = None
    for (c,s,e) in regions:
      R.extend(self.query(c,s,e, extract='raw', **kwargs))
    #efor
    if extract is None:
      return formats.VCF(R)
    #fi

    return formats.VCF.extract(R, extract=extract)
  #edef

  def getVar(self, chromosome, *pargs, **kwargs):
    chromosome = str(chromosome)
    oname = "vcf_%s" % chromosome
    if not self._objectExists(oname):
      utils.error("Could not find chromosome '%s'" % chromosome)
      return None
    #fi
    return self._getObject(oname).getVar(chromosome, *pargs, **kwargs)
  #edef

  def whoHas(self, chromosome, *pargs, **kwargs):
    chromosome = str(chromosome)
    oname = "vcf_%s" % chromosome
    if not self._objectExists(oname):
      utils.error("Could not find chromosome '%s'" % chromosome)
      return None
    #fi
    return self._getObject(oname).whoHas(chromosome, *pargs, **kwargs)
  #edef

  @property
  def phenotypes(self):
    """
    Phenotypic information for the sequenced individuals.
    """
    obj = self._getObject("phen")
    if "famnr" not in obj.columns:
      obj["famnr"] = obj.LLnr.apply(lambda x: str(int(x.split('.')[0])))
    #fi
    return obj
  #edef

  @property
  def percentiles(self):
    """
    Percentile information for all individuals in the LLS study
    """
    return self._getObject("percentiles")
  #edef
  def cgIDToLLNR(self, cgID):
    cgID    = cgID.replace('_240_37-ASM', '')
    possible = self.phenotypes[self.phenotypes.cgID == cgID].LLnr.values 
    if len(possible) > 0:
      return possible[0]
    #fi
    return None
  #edef

  def llnrTpcgID(self, llnr):
    possible = self.self.phenotypes[self.phenotypes.llnr == llnr].cgID.values
    if len(possible) > 0:
      return possible[0]
    #fi
    return None
  #edef

  def passFams(self, pThresh=.90):
    """
    Determine which families pass the novel selection threshold.
    I.e. at least two siblings > pThresh %, and at least one parent > pThresh %
  
    Inputs:
      pThresh: Float. Percentile threshold to satisfy
    Output:
      Set of family IDs (Integers)
    """
    p = self.percentiles
    pF0 = p[(p.percentile >= pThresh) & (p.generation == 1)]
    pF1 = p[(p.percentile >= pThresh) & (p.generation == 2)]
  
    f0Pass = pF0.famnr.unique()
    f1Pass = pF1[ pF1.groupby("famnr").transform(len).llnr >= 2 ].famnr.unique()
  
    famPass = set(f0Pass) & set(f1Pass)
    return famPass
  #edef

  def passIndiv(self, pThresh=.90, relFams=None, sequenced=False):
    """
    Determine, in the families that pass the novel selection threshold, which individuals satisify the criteria.
    Inputs:
      pThresh: Float. Percenfile threshold to satisfy
      relFams: Set of Integers. Family IDs to consider (default is None, then output of passFams is used
      sequenced: Boolean. Only return IDs of individuals that have been sequenced. (Default False
    Output:
      Set of individual IDs.
    """
    # Select all individuals from that family
  
    if relFams is None:
      relFams = self.passFams(pThresh)
    #fi
  
    pPer = self.percentiles
    relIndivs = pPer[pPer.famnr.apply(lambda f: f in relFams) & (pPer.percentile >= pThresh)]
    relIndivs = relIndivs[["famnr", "llnr", "percentile", "motherID", "fatherID", "gender"]]
 
    relIndiv = set(relIndivs.llnr.values)
    if sequenced:
      relIndiv = relIndiv & set(self.phenotypes.LLnr.values)
    #fi
  
    return relIndiv
  #edef


#eclass
