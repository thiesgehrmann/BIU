from .. import utils

np = utils.py.loadExternalModule("numpy")

#################################################################

def order(arr, decreasing=False):
    """
    Return the indexes of the ordering of an array
    """
    o = np.argsort(arr)
    return np.flip(o, axis=0) if decreasing else o
#edef

def pmin(*arrs):
    """
    Parallel minimum
    """
    lens = [ len(arr) for arr in arrs if hasattr(arr, '__len__') ]
    if len(set(lens)) != 1:
        raise Exception
    #fi
    arr_len = lens[0]
    
    arrs = [ arr if hasattr(arr, '__len__') else arr*np.ones(arr_len) for arr in arrs]
    
    par_min = np.array([ min(v) for v in zip(*arrs) ])
    return par_min
#edef

def pmax(maxvalue, arr):
    """
    parallel maximum
    """
    lens = [ len(arr) for arr in arrs if hasattr(arr, '__len__') ]
    if len(set(lens)) != 1:
        raise Exception
    #fi
    arr_len = lens[0]
    
    arrs = [ arr if hasattr(arr, '__len__') else arr*np.ones(arr_len) for arr in arrs]
    
    par_min = np.array([ max(v) for v in zip(*arrs) ])
    return par_min
#edef
#edef


def cummin(arr):
    """
    Return the cululative minimum of an array. (i.e. the next value is always the minimum of the current value and the previous value)
    """
    #ret = reduce(lambda head, nextval: head + [ min(nextval, head[-1])], arr, [max(arr)])[1:]
    ret = arr.copy()
    for i in range(1, len(ret)):
        ret[i] = min(ret[i], ret[i-1])
    #efor
    return np.array(ret)
#edef

def cummax(arr):
    """
    Return the cululative maximum of an array. (i.e. the next value is always the maximum of the current value and the previous value)
    """
    #ret = reduce(lambda head, nextval: head + [ max(nextval, head[-1])], arr, [min(arr)])[1:]
    ret = arr.copy()
    for i in range(1, len(ret)):
        ret[i] = max(ret[i], ret[i-1])
    #efor
    return np.array(ret)
#edef