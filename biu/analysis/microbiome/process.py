from ... import ops
from ... import utils

pd = utils.py.loadExternalModule("pandas")

def relative(A, pseudocount=True):
    """
    Make relative abundances from abundances table
    
    parameters:
    -----------
    
    A: pandas.dataframe
        Columns are 
        Rows are samples
        
    pseudocount: bool|float|int
        If True: add a relative pseudocount
        If False: Do not add a pseudocount
        Otherwise: Add passed as parameter as pseudocount
    """
    if pseudocount is True:
        A = A.add(A.sum(axis=1) / A.sum(axis=1).max(), axis=0) # pseudocount
    elif pseudocount is False:
        pass
    else:
        A = A + float(pseudocount)
    #fi
    RA = A.div(A.sum(axis=1), axis=0)
    
    return RA
#edef

###############################################################################


def group_taxa(A, T, level='genus', exclude=None, exclude_level='taxon_id', A_level='taxon_id'):
    if level == 'species':
        level = 'full_species'
    #fi
    level_order = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'full_species', 'taxon_id']
    level_key   = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']

    if exclude is None:
        exclude = []
    #fi
    if (level not in level_order) or (exclude_level not in level_order) or (A_level not in level_order):
        raise Exception
    #fi
    
    Torig = T.copy()
    Torig['full_species'] = Torig.apply(lambda x: '%s%s%s' % (x.genus, ' ' if x.species is not None else '', x.species if x.species is not None else ''), axis=1)
    T = Torig.copy()

    Aorig = A.copy()
    A = Aorig[[c for c in Aorig.columns if c not in exclude]].copy()
    T = T[T[A_level].isin(A.columns)]
    
    def identifier_at_given_level(t, level):
        idx = level_order.index(level)
        identifier = None
        while identifier is None:
            identifier = t[level_order[idx]]
            idx -= 1
        #ewhile
        return (t[A_level], '%s_%s' % (level_key[idx+1], identifier))
    #edef
    
    group = T.apply(lambda r: identifier_at_given_level(r,level), axis=1)
    group = ops.lst.group(group, key=lambda x: x[1], value=lambda x: x[0])
    puorg = { t: g for (g,T) in group.items() for t in T }

    for (gid, tids) in group.items():
        Aorig[gid] = A[tids].sum(axis=1)
    #efor

    T['taxon_id'] = T[A_level].apply(lambda x: puorg[x])
    T = T.groupby('taxon_id').agg(list)
    for l in level_order:
        if level == l:
            break
        #fi
        T[l] = T[l].apply(lambda x: x[0])
    #efor
    T = T.reset_index()
    T = pd.concat([T, Torig[Torig[A_level].isin(exclude)]])
    
    return Aorig[list(group.keys())+exclude], T.reset_index(drop=True)[level_order]

#edef

###############################################################################

def topn(RA, n=9):
    sel_taxa = RA.mean().sort_values()[-n:].index
    RA_sel = RA[sel_taxa].copy()
    RA_sel['Other'] = RA[[c for c in RA.columns if c not in sel_taxa]].sum(axis=1)
    return RA_sel
#edef
    