from .. import utils

sstats = utils.py.loadExternalModule("scipy.stats")

from collections import namedtuple

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
        oddsratio = 0
        try:
            oddsratio = (t[0][0] / t[0][1]) / (t[1][0] / t[1][1])
        except:
            oddsratio = 100
        #etry
    #fi
    if abcd_values:
        return resTuple(oddsratio, chi2, p, [[a,b],[c,d]], method)
    else:
        return resTuple(oddsratio, chi2, p, table, method)
#edef
