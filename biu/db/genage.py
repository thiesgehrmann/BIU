from ..structures import Dataset2
from .. import utils

pd = utils.py.loadExternalModule('pandas')

class GenAge(Dataset2):
    def __init__(self, *pargs, **kwargs):
        super(GenAge, self).__init__("GenAge", *pargs, **kwargs)
        
        self._obj.add_file('model_genes.csv',
                           utils.Acquire2().curl('https://genomics.senescence.info/genes/models_genes.zip')\
                                     .unzip()\
                                     .select('genage_models.csv'))
        self._obj.add_file('human_genes.csv',
                           utils.Acquire2().curl('https://genomics.senescence.info/genes/human_genes.zip')\
                                     .unzip()
                                     .select('genage_human.csv'))
        
        self._obj.register('model', ['model_genes.csv'],
                           lambda f: pd.read_csv(f['model_genes.csv'],
                                                 index_col=0,
                                                 dtype=dict(zip(['symbol', 'name', 'organism', 'entrez gene id',
                                                                 'avg lifespan change (max obsv)',
                                                                 'lifespan effect', 'longevity influence'],
                                                                 (str,str,str,str,float,str,str)))))
        self._obj.register('human', ['human_genes.csv'],
                           lambda f: pd.read_csv(f['human_genes.csv'],
                                                 index_col=0,
                                                 dtype=dict(zip(['GenAge ID', 'symbol', 'aliases', 'name',
                                                                 'entrez gene id', 'uniprot', 'why', 'band',
                                                                 'location start', 'location end',
                                                                 'orientation', 'acc promoter', 'acc orf',
                                                                 'acc cds', 'references', 'orthologs'],
                                                                 (int,str,str,str,str,str,str,str,int,int,int,
                                                                  str,str,str,str,str)))))
        
    #edef
#eclass