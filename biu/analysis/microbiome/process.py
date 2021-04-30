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
    