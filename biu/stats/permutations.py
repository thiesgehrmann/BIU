from .. import utils

np = utils.py.loadExternalModule('numpy')

def pvalue(statistic, permutation_statistics, n_permutations=None, side='both', pseudocount=1):
    """
    Calculate a permutation test p-value based on a calculated statistic and a set of null-distribution statistics
    
    parameters:
    -----------
    statistic: Float. A calculated statistic
    permutation_statistics. List[Float]|np.array[Float]. A list of null-distribution statistics (e.g. permutations)
    n_permutations: The number of permutations that have been calculated
            (if e.g. you discarded permutations less than the statisic to save memory)
            Default is length of permutation_statistics
    side: 'left' | 'right' | 'both', Calculate exceedences on which side
            left: Count number <= statistic
            right: Count number >= statistic
            both: min(left, right)
    pseudocount: How many pseudocounts to add.
            p-value should never be zero. Adding 1 makes this so.
            
    returns:
    --------
    pvalue: Float
    """
    
    if n_permutations is None:
        n_permutations = len(permutation_statistics)
    #fi
    
    if not isinstance(permutation_statistics, np.ndarray):
        permutation_statistics = np.array(permutation_statistics)
    #fi
    
    exeedences = 0
    
    if side == 'right':
        exeedences = sum(permutation_statistics > statistic)
    elif side == 'left':
        exeedences = sum(permutation_statistics < statistic)
    elif side == 'both':
        exeedences = min(sum(permutation_statistics >= statistic),
                         sum(permutation_statistics <= statistic))
    else:
        raise ValueError("Unknown side: '%s'. See docstring." % side)
    #fi
    
    return (pseudocount + exeedences)/(pseudocount + n_permutations)
#edef