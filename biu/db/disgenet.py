from ..structures import Dataset2
from .. import formats
from .. import utils

class DisGeNet(Dataset2):

    def __init__(self, *pargs, **kwargs):
        """
        Initialize the DisGeNet data structure.
        
        *pargs, **kwargs. See arguments for biu.structures.Dataset2
        
        """
        super(DisGeNet, self).__init__("DisGeNet/", *pargs, **kwargs)
        
        self._obj.add_file('curated_genes.tsv', utils.Acquire2().curl('https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz').gunzip())
        self._obj.add_file('all_genes.tsv', utils.Acquire2().curl('https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_pmid_associations.tsv.gz').gunzip())
        
        self._obj.register('curated_genes', ['curated_genes.tsv'], lambda f: formats.MultiMappingIndex(f['curated_genes.tsv'], sep='\t'))
        self._obj.register('all_genes', ['all_genes.tsv'], lambda f: formats.MultiMappingIndex(f['all_genes.tsv'], sep='\t'))
    #edef
#eclass