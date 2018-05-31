from collections import namedtuple as namedtuple


from ..structures import fileManager as fm
from ..structures import resourceManager as rm

###############################################################################

def urlFileIndex():
  files = {}
  files["geneid2ensemblgene"] = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", "geneid2ensembl.tsv", {})
  files["gene2refseq"]  = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz", "gene2refseq.tsv", {})
  files["geneid2genesymbol"]     = ("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz", "geneinfo.tsv", {})
  files["uniprotmap"]   = ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz", "uniprotmap.tsv", {})
  return { k : (u, 'humanMappings_%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class HumanMapping(fm.FileManager):

  geneid2ensemblgene = None
  #gene2refseq  = None
  geneid2genesymbol     = None
  uniprotmap   = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=["geneid2ensemblgene", "geneid2genesymbol"], **kwargs)

    self.geneid2ensemblgene = rm.TSVMapResourceManager(self, "geneid2ensemblgene", 1, 2, **kwargs)
    self.geneid2genesymbol = rm.TSVMapResourceManager(self, "geneid2genesymbol", 1, 2, **kwargs)

    uniprotFieldNames = ('UniProtKB_AC', 'UniProtKB_ID', 'GeneID', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 
                         'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'PubMed_more')
    #self.ensembltranscript2ensemblgene = rm.TSVMapResourceManager(self, "uniprotmap", 19, 18, fieldNames=uniprotFieldNames, **kwargs)
    #self.ensembltranscript2uniprot     = rm.TSVMapResourceManager(self, "uniprotmap", 19, 1, fieldNames=uniprotFieldNames, **kwargs)
    self.geneid2uniprot = rm.TSVMapResourceManager(self, "uniprotmap", 2, 0, fieldNames=uniprotFieldNames, **kwargs)
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
    geneIDs = self.geneid2ensemblgene.inverse(ensemblID)
    if len(geneIDs) == 0:
      return []
    #fi
    symbols = set.union(*[ set(self.geneid2genesymbol.lookup(g)) for g in geneIDs ])
    return list(symbols)
  #edef

  def getSymbolEnsembl(self, geneSymbol):
    geneIDs = self.geneid2genesymbol.inverse(geneSymbol)
    if len(geneIDs) == 0:
      return []
    #fi
    ensembl = set.union(*[ set(self.geneid2ensemblgene.lookup(gid)) for gid in geneIDs ])
    return list(ensembl)
  #edef

  def getEnsemblGeneID(self, ensemblID):
    geneIDs = self.geneid2ensemblgene.inverse(ensemblID)
    return list(set(geneIDs))
  #edef

  def getGeneIDEnsembl(self, geneID):
    ensemblIDs = self.geneid2ensemblgene.lookup(str(geneID))
    return list(set(ensemblIDs))

  def getSymbolGeneID(self, symbol):
    geneIDs = self.geneid2genesymbol.inverse(symbol)
    return list(set(geneIDs))
  #edef

  def getGeneIDSymbol(self, geneID, other=False):
    geneNameLookup = self.geneid2genesymbol.lookup(str(geneID), withEntry=other)
    if other:
      geneNames = [ r[0] for r in geneNameLookup ]
      otherEntries = [ n for r in geneNameLookup for n in self.geneid2ensemblgene[r[1]][5].split('|') if n != '-' ]
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
    ensemblID = None if len(ensemblID) == 0 else symbol[0]

    return self.__allIDs(geneID, ensemblID, symbol)
  #edef

#eclass
