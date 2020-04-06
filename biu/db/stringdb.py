from ..structures import Dataset2
from .. import utils

np = utils.py.loadExternalModule('numpy')
pd = utils.py.loadExternalModule('pandas')

class StringDB(Dataset2):
    """
    The String DataBase
    Retrieve the network around a specified set of genes.
    
    Example usage:
    --------------
    string = biu.db.StringDB(9606) # Human
    
    net = string.get(['FOXO3', 'HIF1A', 'TNF', 'SOD2'])
    """
    def __init__(self, organism=9606, *pargs, **kwargs):
        """
        Initialize a string DB instance.
        
        parameters:
        -----------
        organism: int
            The Taxonomic ID of the organism you want to use
            9606 : Homo sapiens
        """
        super(StringDB, self).__init__("StringDB/%d" % organism, *pargs, **kwargs)
        
        object.__setattr__(self, 'organism', int(organism))
        object.__setattr__(self, '_thresholds', {})
        
        url_head = 'https://stringdb-static.org/download'
        self._obj.add_file('protein_info',
                           utils.Acquire2()\
                               .curl('%s/protein.info.v11.0/%d.protein.info.v11.0.txt.gz' % (url_head, self.organism))\
                               .gunzip()\
                               .cmd("sed {} -e 's/^%d[.]//' -e 's/ %d[.]/ /'" % (self.organism, self.organism), '{}'))
        self._obj.add_file('protein_links',
                           utils.Acquire2()\
                               .curl('%s/protein.links.detailed.v11.0/%d.protein.links.detailed.v11.0.txt.gz' % (url_head, self.organism))\
                               .gunzip()\
                               .cmd("sed {} -e 's/^%d[.]//' -e 's/ %d[.]/ /'" % (self.organism, self.organism), '{}'))
        
        self._obj.register('_edges', ['protein_links'], lambda f: pd.read_csv(f['protein_links'], sep=' '))
        self._obj.register('_nodes', ['protein_info'], lambda f: pd.read_csv(f['protein_info'], sep='\t'))
    #edef
    
    def get(self, selection, extend=0, threshold=400, self_links=True):
        """
        Get the network for a specified set of genes.

        parameters:
        -----------

        selection: List[String]
            List of Ensembl protein IDs of HGNC symbols you want to consider
        extend: Integer
            Extend the selection with the n'th neighbour (default=0)
        threshold: Integer [0-1000]
            Select only connections between genes that have an edge in any network stronger than this threshold
        self_links: Boolean
            Include self-links in the network
            These don't exist in the original STRING db, but we add them to retain the
            definition of the original selection. E.g. if there are no links associated with the gene, then we 

        threshold: float [0-1000]
            The threshold to use on the network
            For the human, it changes the number of links in this way:
              0   -> 11759454
              100 -> 11759454
              200 ->  6894162
              300 ->  3256482
              400 ->  2000556
              600 ->  1045444
              700 ->   841068
              800 ->   728090
              900 ->   648304

        Returns: pd.DataFrame

        """

        selection = set(selection)
        all_genes = set(self._nodes.protein_external_id) | set(self._nodes.preferred_name)
        good_selection = selection & all_genes
        if len(good_selection) != len(selection):
            biu.utils.msg.warning('Could not find the following genes:\n%s' % ', '.join(selection - all_genes) )
        #fi

        sel_ids = self._nodes[self._nodes.protein_external_id.isin(good_selection) | \
                              self._nodes.preferred_name.isin(good_selection)].protein_external_id

        orig_sel_ids = sel_ids


        if threshold not in self._thresholds:
            self._thresholds[threshold] = (self._edges[self._edges.columns[2:]] >= threshold).any(axis=1)
        #fi

        threshnet = self._edges[self._thresholds[threshold]]

        while extend > 0:
            print(extend)
            net_sel = threshnet[threshnet.protein1.isin(sel_ids) | threshnet.protein2.isin(sel_ids)]
            sel_ids = set(net_sel.protein1) | set(net_sel.protein2)
            extend = extend - 1
        #ewhile

        net_sel = threshnet[threshnet.protein1.isin(sel_ids) & threshnet.protein2.isin(sel_ids)]

        if self_links:
            self_links = [ (p, p, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000) for p in orig_sel_ids ]
            self_links = pd.DataFrame(self_links, columns=net_sel.columns)
            net_sel = pd.concat([net_sel, self_links])
        #fi

        return net_sel
    #edef  
#eclass