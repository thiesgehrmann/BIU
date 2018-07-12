from collections import namedtuple as namedtuple

from ..config import settings
from ..structures import Dataset
from .. import formats
from .. import utils

###############################################################################

class HumanMapping(Dataset):

  geneid2ensemblgene = None
  #gene2refseq  = None
  geneid2genesymbol     = None
  uniprotmap   = None

  def __init__(self, **kwargs):
    fileIndex = self.__genFileIndex()
    Dataset.__init__(self, fileIndex)

    self._registerObject("geneid2ensemblgene", formats.TSVMap, ["geneid2ensemblgene"], fileIndex["geneid2ensemblgene"].path, 1, 2, delimiter='\t')
    self._registerObject("geneid2genesymbol", formats.TSVMap, ["geneid2genesymbol"], fileIndex["geneid2genesymbol"].path, 1, 2, delimiter='\t')

    uniprotFieldNames = ('UniProtKB_AC', 'UniProtKB_ID', 'GeneID', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR',
                         'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'PubMed_more')

    self._registerObject("geneid2uniprot", formats.TSVMap, ["uniprotmap"], fileIndex["uniprotmap"].path, 1, 2, names=uniprotFieldNames)

  #edef

  def __genFileIndex(self, where=None):
    files = {}
    finalPath = '%s/maps/humanMappings' % (settings.getDataDir() if where is None else where)

    files["geneid2ensemblgene"] = utils.Acquire(where=where).curl("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz").gunzip().finalize('%s/geneid2ensembl.tsv' % finalPath)
    files["gene2refseq"]        = utils.Acquire(where=where).curl("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz").gunzip().finalize('%s/gene2refseq.tsv' % finalPath)
    files["geneid2genesymbol"]  = utils.Acquire(where=where).curl("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz").gunzip().finalize('%s/geneinfo.tsv' % finalPath)
    files["uniprotmap"]         = utils.Acquire(where=where).curl("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz").gunzip().finalize('%s/uniprotmap.tsv' % finalPath)
    return files
  #edef

  ############################################################################### 

  @property
  def __geneid2ensemblgene(self):
    return self._getObject("geneid2ensemblgene")
  #edef

  @property
  def __geneid2genesymbol(self):
    return self._getObject("geneid2genesymbol")
  #edef

  ###############################################################################

  def getEnsemblUniprot(self, ensemblID):
    geneIDs = self.getEnsemblGeneID(ensemblID)
    if len(geneIDs) == 0:
      return []
    #fi
    uniprot = set.union(*[ set(self.geneid2uniprot.lookup(g)) for g in geneIDs ])
    return uniprot
  #edef

  def getUniprotEnsembl(self, uniprot):
    geneIDs = self.geneid2uniprot.inverse(uniprot)
    if len(geneIDs) == 0:
      return []
    #fi
    ensembl = set.union(*[ set(self.geneid2ensembl.lookup(g)) for g in geneIDs ])
    return ensembl
  #edef
    

  def getEnsemblSymbol(self, ensemblID):
    geneIDs = self.__geneid2ensemblgene.inverse(ensemblID)
    if len(geneIDs) == 0:
      return []
    #fi
    symbols = set.union(*[ set(self.__geneid2genesymbol.lookup(g)) for g in geneIDs ])
    return list(symbols)
  #edef

  def getSymbolEnsembl(self, geneSymbol):
    geneIDs = self.__geneid2genesymbol.inverse(geneSymbol)
    if len(geneIDs) == 0:
      return []
    #fi
    ensembl = set.union(*[ set(self.__geneid2ensemblgene.lookup(gid)) for gid in geneIDs ])
    return list(ensembl)
  #edef

  def getEnsemblGeneID(self, ensemblID):
    geneIDs = self.__geneid2ensemblgene.inverse(ensemblID)
    return list(set(geneIDs))
  #edef

  def getGeneIDEnsembl(self, geneID):
    ensemblIDs = self.__geneid2ensemblgene.lookup(str(geneID))
    return list(set(ensemblIDs))

  def getSymbolGeneID(self, symbol):
    geneIDs = self.__geneid2genesymbol.inverse(symbol)
    return list(set(geneIDs))
  #edef

  def getGeneIDSymbol(self, geneID, other=False):
    geneNameLookup = self.__geneid2genesymbol.lookup(str(geneID), withEntry=other)
    if other:
      geneNames = [ r[0] for r in geneNameLookup ]
      otherEntries = [ n for r in geneNameLookup for n in self.__geneid2ensemblgene[r[1]][5].split('|') if n != '-' ]
      return list(set(geneNames + otherEntries))
    else:
      return list(set(geneNameLookup))
    #fi
  #edef

  __allIDs = namedtuple("GeneIDMapping", ['geneID', 'ensemblID', 'symbol'])

  def fromEnsembl(self, ensemblID):
    geneID = self.getEnsemblGeneID(ensemblID)
    geneID = None if len(geneID) == 0 else geneID[0]
    symbol = self.getEnsemblSymbol(ensemblID)
    symbol = None if len(symbol) == 0 else symbol[0]

    return self.__allIDs(geneID, ensemblID, symbol)
  #edef

  def fromGeneID(self, geneID):
    ensemblID = self.getGeneIDEnsembl(geneID)
    ensemblID = None if len(ensemblID) == 0 else ensemblID[0]
    symbol = self.getGeneIDSymbol(geneID)
    symbol = None if len(symbol) == 0 else symbol[0]

    return self.__allIDs(geneID, ensemblID, symbol)
  #edef

  def fromSymbol(self, symbol):
    geneID = self.getSymbolGeneID(symbol)
    geneID = None if len(geneID) == 0 else geneID[0]
    ensemblID = self.getSymbolEnsembl(symbol)
    ensemblID = None if len(ensemblID) == 0 else ensemblID[0]

    return self.__allIDs(geneID, ensemblID, symbol)
  #edef

#eclass
