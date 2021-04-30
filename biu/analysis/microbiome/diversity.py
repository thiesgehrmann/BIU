from ... import utils
from ... import ops

np  = utils.py.loadExternalModule('numpy')

###############################################################################

def shannon_entropy(D, assume_lowest_value_equals_zero=True):
    """
    Calculate the shannon entropy
    
    parameters:
    -----------
    D: array of floats
    assume_lowest_value_equals_zero : boolean
        We assume that the lowest value is actually a zero, and is just a pseudocount
        These values are ignored in the calculation
    """
    if assume_lowest_value_equals_zero:
        m = min(D)
        D = [ d for d in D if d > m ]
    #fi
    return - sum(np.multiply(D, np.log(D)))
#edef

###############################################################################

def partial_shannon_entropy(D, pos=None, assume_lowest_value_equals_zero=True):
    """
    Calculate the shannon entropy
    
    parameters:
    -----------
    D: array of floats
    pos: int
        Which position of the diversity to use? STARTING AT THE MOST ABUNDANT
    assume_lowest_value_equals_zero : boolean
        We assume that the lowest value is actually a zero, and is just a pseudocount
        These values are ignored in the calculation
    """
    if assume_lowest_value_equals_zero:
        m = min(D)
        D = [ d for d in D if d > m ]
    #fi
    
    D = np.sort(D)
    norm = np.cumsum(D)
    E = np.zeros_like(D)
    for i, n in enumerate(norm):
        E[i] = shannon_entropy(D[:i+1] / norm[i], not assume_lowest_value_equals_zero)
    #efor
    
    
    if pos is None:
        return E
    elif pos < len(E):
        return E[pos]
    else:
        return 0
    #fi
#edef

###############################################################################

def area_under_partial_shannon_entropy(D):
    E = partial_shannon_entropy(D)
    nz = sum(D == min(D))
    return sum(E[nz:]) / (len(D)-nz)
#edef

###############################################################################


def kl_divergence_discrete(P,Q):
    return - np.sum(np.multiply(P, np.log2(np.divide(P,Q))))
#edef

# For 1000s of samples, this is too slow!!!!
# Maybe I just need to do it!
def jensen_shannon_divergence(A,B):
    return (kl_divergence_discrete(A,B) + kl_divergence_discrete(B,A))/2
#edef

def jensen_shannon_divergence(A,B):
    return shannon_entropy((A+B)/2) - (shannon_entropy(A) + shannon_entropy(B))/2
#edef

def jensen_shannon_distance(A,B):
    return np.sqrt(jensen_shannon_divergence(A,B))
#edef