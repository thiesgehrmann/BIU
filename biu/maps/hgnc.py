from ..structures import Dataset
from .. import formats
from .. import utils
from ..config import settings

###############################################################################


class HGNC(Dataset):

  versions = {}

  def __init__(self, where=None):

    fileIndex = self.__genFileIndex()
    Dataset.__init__(fileIndex)

    indexes = [ ('hgncID', 0),
                ('symbol', 1),
                ('entrezID', 18),
                ('ensemblGeneID' : 19),
                ('ucscID', 20)]

    for (idxName, idx) in indexes:
      self._registerObject("_map_%s" % idxName, formats.TSVIndex, fileIndex["hgnc_db"], idx)
    #efor
  #edef

  def __genFileIndex(self, where=None):
    finalPath = '%s/maps/hgnc' % ( (settings.getDataDir() if where is None else where), version)
    files = {}
    files["hgnc_db"]  = utils.Acquire().curl("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt").finalize('%s/hgnc_complete_set.tsv' % finalPath)
    return files
  #edef

  ###############################################################################

  def fromHGNC(self, idx):
    return self._getObject('_map_hgncID')[idx]
  #edef

  def fromSymbol(self, idx):
    return self._getObject('_map_symbol')[idx]
  #edef

  def fromEnsembl(self, idx):
    return self._getObject('_map_ensemblGeneID')[idx]
  #edef

  def fromEntrez(self, idx):
    return self._getObject('_map_entrezID')[idx]
  #edef

  def from
