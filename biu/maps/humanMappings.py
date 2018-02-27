
from ..structures import fileManager as fm
from ..structures import resourceManager as rm

###############################################################################

def urlFileIndex():
  files = {}
  files["gene2ensembl"] = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", "gene2ensembl.tsv", {})
  files["gene2refseq"]  = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz", "gene2refseq.tsv", {})
  files["geneInfo"]     = ("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz", "geneinfo.tsv", {})
  files["uniprotmap"]   = ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz", "uniprotmap.tsv", {})
  return { k : (u, 'humanMappings_%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class HumanMapping(fm.FileManager):

  gene2ensembl = None
  #gene2refseq  = None
  geneInfo     = None
  uniprotmap   = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=["gene2ensembl", "geneInfo"], **kwargs)

    self.gene2ensembl = rm.TSVMapResourceManager(self, "gene2ensembl", 1, 2, **kwargs)
    self.geneInfo = rm.TSVMapResourceManager(self, "geneInfo", 1, 2, **kwargs)
  #edef

  ############################################################################### 

  def getEnsemblSymbol(self, ensemblID):
    geneIDs = self.gene2ensembl.inverse(ensemblID)
    if len(geneIDs) == 0:
      return []
    #fi
    symbols = set.union(*[ set(self.geneInfo.lookup(g)) for g in geneIDs ])
    return list(symbols)
  #edef

  def getSymbolEnsembl(self, geneSymbol):
    geneIDs = self.geneInfo.inverse(geneSymbol)
    if len(geneIDs) == 0:
      return []
    #fi
    ensembl = set.union(*[ set(self.gene2ensembl.lookup(gid)) for gid in geneIDs ])
    return list(ensembl)
  #edef

  def getEnsemblGeneID(self, ensemblID):
    geneIDs = self.gene2ensembl.inverse(ensemblID)
    return list(set(geneIDs))
  #edef

  def getGeneIDEnsembl(self, geneID):
    ensemblIDs = self.gene2ensembl.lookup(geneID)
    return list(set(ensemblIDs))

  def getSymbolGeneID(self, symbol):
    geneIDs = self.geneInfo.inverse(symbol)
    return list(set(geneIDs))
  #edef

  def getGeneIDSymbol(self, geneID, other=False):
    geneNameLookup = self.geneInfo.lookup(geneID, withEntry=other)
    if other:
      geneNames = [ r[0] for r in geneNameLookup ]
      otherEntries = [ n for r in geneNameLookup for n in self.gene2ensembl[r[1]][5].split('|') if n != '-' ]
      return list(set(geneNames + otherEntries))
    else:
      return list(set(geneNameLookup))
    #fi
  #edef

#eclass
