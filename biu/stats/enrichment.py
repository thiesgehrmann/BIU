from .. import utils
from .. import ops

sstats = utils.py.loadExternalModule("scipy.stats")
np     = utils.py.loadExternalModule('numpy')

plt = utils.py.loadExternalModule('matplotlib.pylab') 
fc  = utils.py.loadExternalModule('fastcluster')
sp  = utils.py.loadExternalModule('scipy')

from collections import namedtuple

from . import permutations

@utils.decorators.deprecated("setEnrichment is deprecated. Use set_enrichment instead.")
def setEnrichment(your_set, other_set, universe):
    """
    Perform set enrichment using either a fisher exact test or the chi2 test.
    parameters:
    -----------
    your_set:  list. Elements you want to test for enrichment
    other_set: list. Elements you want to see whether they are enriched in your_set
    universe:  list. Total universe of elements
    abcd_values: Boolean. If True, it will return the actual element values in the contingency table, rather than just counts
    
    returns:
    Named tuple with:
        * oddsratio: fisher oddsratio
        * c2statistic : chi2 test statistic
        * pvalue : pvalue of test
        * table: contingency table [ [a,b],[c,d] ]
           - a: Overlap of the two sets
           - b: What is in other_set but not in your_set
           - c: what is in your_set but not in other_set
           - d: What is in universe but not in your_set or other_set
        * table_values: contingency table values [ [a,b],[c,d] ]
           - As see above
        * method : fisher|c2statistic
    """
    return set_enrichment(your_set, other_set, universe)
#edef

def set_enrichment(your_set, other_set, universe):
    """
    Perform set enrichment using either a fisher exact test or the chi2 test.
    parameters:
    -----------
    your_set:  list. Elements you want to test for enrichment
    other_set: list. Elements you want to see whether they are enriched in your_set
    universe:  list. Total universe of elements
    abcd_values: Boolean. If True, it will return the actual element values in the contingency table, rather than just counts
    
    returns:
    Named tuple with:
        * oddsratio: fisher oddsratio
        * c2statistic : chi2 test statistic
        * pvalue : pvalue of test
        * table: contingency table [ [a,b],[c,d] ]
           - a: Overlap of the two sets
           - b: What is in other_set but not in your_set
           - c: what is in your_set but not in other_set
           - d: What is in universe but not in your_set or other_set
        * table_values: contingency table values [ [a,b],[c,d] ]
           - As see above
        * method : fisher|chi2
    """

    
    resTuple = namedtuple("setEnrichmentResult", [ 'oddsratio', 'c2statistic', 'pvalue', 'table', 'table_values', 'method'])

    universe  = set(universe)
    your_set  = set(your_set) & universe
    other_set = set(other_set) & universe
    
    a = your_set & other_set
    b = other_set - your_set
    c = your_set - other_set
    d = universe - (your_set | other_set)
    
    table = [ [len(a), len(b)], [len(c), len(d)]]
    if min(min(table)) <= 5:
        method = 'fisher'
        oddsratio, p = sstats.fisher_exact(table)
        chi2 = None
    else:
        method = 'chi2'
        chi2, p, dof, expected = sstats.chi2_contingency(table)
        oddsratio = 100
        if table[1][0] > 0 and table[0][1] > 0:
            oddsratio = table[0][0] * table[1][1] / (table[1][0] * table[0][1])
        else:
            oddsratio = np.inf
        #fi
    #fi

    return resTuple(oddsratio, chi2, p, table, [[a,b],[c,d]], method)
#edef

def gsea(scores, membership, sort=True, sort_abs=True, p=1, side='both',
         max_perm=1000, min_perm=100, perm_thresh=0.2, plot=None):
    """
    Gene Set Enrichment Analysis.
    
    parameters:
    -----------
    scores:      A list of scores. Each score refers to a gene/locus/object/whatever.
                 NOTE: the SMALLEST score will be at the TOP of the list. Thus, a ranking of
                                    [ 4,2,1,5 ] -> [ 1,2,4,5 ]
    membership:  A list of 0/1 for each object, indicating whether it is in the desired set or not.
    sort:        Are the items already sorted?
    sort_abs:    Boolean. If sort, then sort the scores with absolute value (or not)
    p:           Float, An exponent p to control the weight of the step.
    side:        'left' | 'right' | 'both', Calculate exceedences on which side
                   left: Count number <= statistic
                   right: Count number >= statistic
                   both: min(left, right)
    max_perm:    Integer. The maximum number of permutations to perform when calculating p-values
                 We attempt to prevent many unnecessary permutations.
    min_perm:    Integer. The absolute minimum number of permutations to perform.
    perm_thresh: Float. Only perform more permutations than min_perm if the % of exceedences of the
                        statistic is less than this value
    plot:        None|matplotlib.axis.
                 if not None, plot the histogram 
    
    Returns:
    --------
    a named tuple with:
    (es=enrichment_score,
     p=pvalue,
     i=number_of_genes_at_peak,
     idx=original_index_of_set_genes_at_peak)
    """
    
    nt = namedtuple('GSEA_Result', [ 'es', 'p', 'i', 'idx'])
    
    if len(scores) != len(membership):
        raise ValueError("gsea: scores and membership must be same length")
    #fi
    
    L = zip(range(len(scores)), scores, membership)
    L = sorted(L, key=lambda o: o[1])
    
    O, S, M = zip(*L)
    S = np.array(S)
    M = 1*np.array(M)
    O = np.array(O)

    N  = len(S)
    NH = sum(M)
    
    if NH == 0:
        return nt(0, 1.0, N, [])
    #fi
    
    def calculate_es(s, m):
        NR = sum(s[m==1]**p)

        def pmiss(i):
            return sum(m[:i]==0) / (N - NH)
        #edef
        def phit(i):
            return sum((s[:i][m[:i]==1])**p)/NR
        #edef
        
        es_i = [ (i, phit(i) - pmiss(i)) for i in range(1, N) ]
        
        i, es = sorted(es_i, key=lambda x: x[1])[-1]
        
        return es, es_i, i
    #edef
    
    es, es_i, i = calculate_es(S, M)
    index_i = O[np.where(M[:i] == 1)]
    
    perm_es = [ calculate_es(S, np.random.choice(M, N, replace=False))[0] for i in range(min_perm) ]
    perm_steps = int(np.ceil(max_perm / 10))
    
    nex = 0
    while len(perm_es) < max_perm:
        if side == 'left':
            nex = len([e for e in perm_es if e <= es ])
        elif side == 'right':
            nex = len([e for e in perm_es if e >= es ])
        elif side == 'both':
            nex = min(len([e for e in perm_es if e <= es ]), 
                      len([e for e in perm_es if e >= es ]))
        else:
            raise ValueError("Unknown side: '%s'. See docstring." % side)
        #fi
        
        if nex / len(perm_es) > perm_thresh:
            break
        #fi
        
        perm_es.extend([ calculate_es(S, np.random.choice(M, N, replace=False))[0] for i in range(perm_steps) ])
        
    #ewhile
    return nt(es, permutations.pvalue(es, perm_es, side=side), i, index_i)
#edef

##############################################################################

class EnrichmentNetwork(object):
    """
    Make a network visualization from the results of an enrichment.
    """
    def __init__(self, enrichments, q_col='q', table_col='table', table_values_col='table_values'):
        """
        
        parameters:
        -----------
        enrichments: pandas.DataFame
          Output from an enrichment. e.g. biu.db.KEGG.enrich, or biu.db.Reactome.enrich.
        q_col: String
            The name of the column with the corrected p-values
        table_col: String
            The name of the column with the contingency table
        table_values_col: String
            The name of the column with the object names in  the contingency table (gene names rather than counts)
          
        Properties:
        -----------
        nodes: the original enrichments
        edges: The distance between terms
        
        draw(): Draw the network.
        """
        def distance(a,b):
            if a == b:
                return 0
            #fi
            return 1/(np.log10(len(a&b) + 1)+1)
        #edef

        feat = enrichments[[q_col, table_col, table_values_col]].rename(columns={
            q_col: 'q',
            table_col : 'table',
            table_values_col : 'table_values'}).copy()
        feat['total'] = feat.table_values.apply(lambda x: x[0][0] | x[0][1])
        feat['n'] = feat.table.apply(lambda x: x[0][0])

        values = feat.total.to_dict()

        D = np.zeros([feat.shape[0]]*2)

        for i, (node_i, values_i) in enumerate(values.items()):
            for j, (node_j, values_j) in enumerate(values.items()):
                D[i,j] = D[j,i] = distance(values_i, values_j)
            #efor
        #efor
        
        self._enrichments = enrichments
        self._nodes = feat
        self._edges = D
    #edef
    
    @property
    def nodes(self):
        return self._nodes
    #edef
    
    @property
    def edges(self):
        return self._edges
    #edef
    
    def _embed(self, network):
        from sklearn.manifold import MDS
        E = MDS(n_components=2, dissimilarity='precomputed').fit_transform(self.edges)
        return E
    #edef
    
    def draw(self, distance_threshold=0.3, ax=None, cmap=plt.get_cmap('plasma'), nodes=None, n_clusters=1, min_qval=None):

        from sklearn.manifold import MDS
        from sklearn.manifold import Isomap
        import fastcluster
        import scipy as sp

        E = self._embed(self._edges)

        self._nodes['x'] = E[:,0]
        self._nodes['y'] = E[:,1]

        if ax is None:
            fig, axes = utils.figure.subplots(figsize=(10,10), dpi=300)
            ax = axes[0]
        #fi
        
        if min_qval is None:
            min_qval = self._nodes.q.min()
        #fi

        color = self._nodes.q.apply(lambda x: -np.log10(x)) / -np.log10(min_qval)
        ax.scatter(self._nodes.x, self._nodes.y, s=self._nodes.n*4, edgecolor='k', c=cmap(color))
        for i, r in self._nodes.iterrows():
            ax.text(r.x, r.y, str(i))
        #efor

        #L = fc.linkage(sp.spatial.distance.squareform(self.edges, checks=False), method='complete')
        L = fc.linkage(E, metric='euclidean', method='complete')
        clusters = ops.lst.flatten(sp.cluster.hierarchy.cut_tree(L, n_clusters=n_clusters))

        for i,j in np.ndindex(self.edges.shape):
            if i > j:
                continue
            #fi

            if self.edges[i,j] < distance_threshold:
                if clusters[i] == clusters[j]:
                     ax.plot([E[i,0],E[j,0]], [E[i,1],E[j,1]], zorder=-1, c='#565656', alpha=0.5)
                else:
                    ax.plot([E[i,0],E[j,0]], [E[i,1],E[j,1]], zorder=-2, c='#eaeaea')
                #fi
            #fi
        #efor
        
        xlim, ylim = (ax.get_xlim(), ax.get_ylim())
        plotlim = ylim[0] + (ylim[1]-ylim[0])/10
        
        sizes = np.array([ 10, 25, 50, 100, 200, 300 ])
        sizes = sizes[ sizes <= max(self._nodes.n) ]
        stepsize = (ylim[1]-ylim[0])/10 / len(sizes)
        
        plt.scatter([xlim[0]+(xlim[1]-xlim[1])*0.1] * len(sizes),
                    ylim[0] + (1 + np.array(range(len(sizes))))[::-1]*stepsize,
                    c='k', edgecolor='k', s=sizes*4)
        
        dendfig, dendaxes = utils.figure.subplots(figsize=(20,5))
        dend = sp.cluster.hierarchy.dendrogram(L, ax=dendaxes[0], labels=self.nodes.index, )

        return ax.get_figure()
    #edef
#eclass
