from .. import utils

sstats = utils.py.loadExternalModule("scipy.stats")

from collections import namedtuple

def setEnrichment(yourSet, otherSet, allAnnotatedObjects):

    resTuple = namedtuple("setEnrichmentResult", [ 'oddsratio', 'c2statistic', 'pvalue', 'table', 'method' ])

    allAnnotatedObjects = set(allAnnotatedObjects)
    yourSet  = set(yourSet) & allAnnotatedObjects
    otherSet = set(otherSet) & allAnnotatedObjects
    
    a = yourSet & otherSet
    b = otherSet - yourSet
    c = yourSet - otherSet
    d = allAnnotatedObjects - (yourSet | otherSet)
    
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
    return resTuple(oddsratio, chi2, p, table, method)
#edef
