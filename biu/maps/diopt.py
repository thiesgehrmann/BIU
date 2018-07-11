from ..structures import Dataset
from .. import formats
from ..config import settings
from .. import utils

from collections import namedtuple

pd = utils.py.loadExternalModule("pandas")

class DIOPT(Dataset):
  
  __speciesIDS = {'h_sapiens' : 9606,
        's_pombe': 4896,
        's_cerevisiae' : 4932,
        'c_elegans' : 6239,
        'd_melanogaster' : 7227,
        'd_reiro' : 7955,
        'x_tropicalis' : 8364,
        'm_musculus' : 10090,
        'r_norvegicus' : 10116 }

  __speciesNames = {'h_sapiens' : 'homo sapiens',
        's_pombe': 'schizosaccharomyces pombe',
        's_cerevisiae' : 'saccharomyces cerevisiae',
        'c_elegans' : 'caenorhabditis elegans',
        'd_melanogaster' : 'drosophila melanogaster',
        'd_reiro' : 'danio rerio',
        'x_tropicalis' : 'xenopus tropicalis',
        'm_musculus' : 'mus musculus',
        'r_norvegicus' : 'rattus norvegicus' }
  
  species = list(__speciesIDS.keys())
  
  __fields = [ 'species', 'ncbi', 'symbol', 'id', 'score', 'confidence', 'database']
  __namedtupleObject = namedtuple('DIOPTResult', __fields)
  
  def __init__(self, genome=list(__speciesIDS.keys())[0], where=None, **kwargs):
    finalPath = '%s/diopt' % (settings.getWhere() if where is None else where)
    files = {}
    files['data'] = utils.Acquire(where=where).touch().finalize('%s/%s.data.sqlite' % (finalPath, genome))
    
    Dataset.__init__(self, files, **kwargs)
    self._registerObject('_data', formats.SQLDict, [ 'data' ], files['data'].path)
    
    self.__genome = self.__speciesIDS[genome]
  #edef
  
  def batch(self, queryIDs, info=['symbol', 'id']):
    D = []
    for q in queryIDs:
      res  = self[q]
      resRow = []
      for s in self.species:
        rel = [ r for r in res if (r.species == s) and (r.confidence != 'low') ]
        if len(rel) > 0:
          resRow.append(tuple( [ getattr(rel[0],f) for f in info] ))
        #fi
        else:
          resRow.append(tuple([]))
        #fi
      #efor
      D.append(resRow)
    #efor
    return pd.DataFrame(D, columns=self.species, index=queryIDs)
  #edef
  
    
  
  def __getitem__(self, queryID):
    return self.__lookup(queryID)
  #edef
  
  def __delitem__(self, queryID):
    del self._data[str(queryID)]
  #edef
  
  def __lookup(self, queryID):
    queryID = str(queryID)
    if queryID not in self._data:
      self._data[queryID] = self.__query(queryID)
    #fi
    r = self._data[queryID]
    return [ self.__namedtupleObject(*r) for r in self._data[queryID] ]
  #edef
    
  
  def __query(self, queryID):
    ao = utils.Acquire().curl('http://www.gene2function.org/search/get_ortholog/%d/%s/best_match' % (self.__genome, queryID))
    file = ao.acquire()
    return self.__parseQuery('\n'.join(open(file, 'r').readlines()))
  #edef
  
  def __parseQuery(self, queryResult):
    rowFields = ['ncbi', 'symbol', 'hdc', 'species', 'id', 'database', 'diopt_score', 'best_score',
          'best_score_reverse', 'confidence']
    nFields = len(rowFields)
    result = []
    for row in utils.html.table2csv(queryResult):
      row = row[:nFields]
      if len(row) < nFields:
        continue
      #fi
      row = [ f.strip() for f in row ]
      row = dict(zip(rowFields, row))
      thisSpecies = [ name for name in self.__speciesNames if self.__speciesNames[name] in row['species'].lower()]
      if len(thisSpecies) != 1:
        continue
      #fi
      row['species'] = thisSpecies[0]
      row['score'] = row['diopt_score']
      result.append([ row[k] for k in self.__fields ])
    #efor
    return result
  #edef
#eclass
