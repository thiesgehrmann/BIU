import scipy.stats as sstats

def setEnrichment(yourSet, otherSet, allAnnotatedObjects):
    allAnnotatedObjects = set(allAnnotatedObjects)
    yourSet  = set(yourSet) & allAnnotatedObjects
    otherSet = set(otherSet) & allAnnotatedObjects
    
    a = yourSet & otherSet
    b = otherSet - yourSet
    c = yourSet - otherSet
    d = allAnnotatedObjects - (yourSet | otherSet)
    
    table = [ [len(a), len(b)], [len(c), len(d)]]
    if min(min(table)) <= 5:
        return sstats.fisher_exact(table), table, 'fisher'
    else:
        return sstats.chi2_contingency(table), table, 'chi2'
    #fi
#edef
