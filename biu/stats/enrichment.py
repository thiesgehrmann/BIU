from .. import utils

sstats = utils.py.loadExternalModule("scipy.stats")
np     = utils.py.loadExternalModule('np')

from collections import namedtuple

from . import permutations

@utils.decorators.deprecated("setEnrichment is deprecated. Use set_enrichment instead.")
def setEnrichment(your_set, other_set, universe, abcd_values=False):
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
        * method : fisher|c2statistic
    """
    return set_enrichment(your_set, other_set, universe, abcd_values)
#edef

def set_enrichment(your_set, other_set, universe, abcd_values=False):
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
        * method : fisher|chi2
    """

    
    resTuple = namedtuple("setEnrichmentResult", [ 'oddsratio', 'c2statistic', 'pvalue', 'table', 'method'])

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
    if abcd_values:
        return resTuple(oddsratio, chi2, p, [[a,b],[c,d]], method)
    else:
        return resTuple(oddsratio, chi2, p, table, method)
    #fi
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