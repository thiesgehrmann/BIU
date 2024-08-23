import hashlib
import math

from .. import utils

np = utils.py.loadExternalModule("numpy")

def freq(lst, key=lambda x: x):
  """
    Count the number of occurances in a list.
    Input:
        lst: an iterable list
        key: a callable function that is called upon the value, the return value of which is counted
    
    Output: Dictionary with count of each item
  """
  F = {}
  for item in lst:
    F[key(item)] = F.get(key(item), 0) + 1
  #efor
  return F
#edef

def qerf(lst, key=lambda x: x):
    """
    Return the frequency of objects in a list, indexed by their frequency, rather than by the item itself.
    parameters:
    -----------
    lst: list[obj]
        List of hashable objects
    key: callable(obj)
        A function which will be called upon the objects in the list
        
    example:
    qerf([1,2,3,3,4,4,4])
     { 1: [1,2],
       2: [3],
       3: [4]
       }
    """
    f = freq(lst, key)
    q = { }
    for (i,c) in f.items():
        q[c] = q.get(c,[]) + [i]
    #efor
    return q
#edef

def group(lst, key=lambda x: x[0], value=lambda x: x):
  """
    Group items based on a certain key.
    Input: lst: an iterable list of tuple (or indexable values)
           key: A function to determine grouping
    Output: Dictionary with items grouped by key.
  """
  F = {}
  for item in lst:
    k = key(item)
    F[k] = F.get(k, []) + [value(item)]
  #efor
  return F
#edef

def flatten(lst):
  """
  flatten: flatten a list of lists of elements into a list of elements
  Input: A list of lists of elements
  Output: A list of elements
  """
  return [ item for group in lst for item in group ]
#edef

def hash(arr, strategy="tmb", f=hashlib.md5):
  """
  hash: Hash a list of elements
  Input: arr: A list of elements with a __str__ attribute
         strategy: tmb (sparse selection. selects 10 elements from list at top, middle and bottom for hash)
                   all (Full selection)
         f: Hash function (Default is hashlib.md5)
  output: Hex digest of hash
  """
  h = f()

  if strategy == 'tmb': # Top Middle Bottom
    middle = int(len(arr) / 2)
    for o in arr[:10] + arr[middle:middle+10] + arr[:-10]:
      h.update(str(o).encode())
    #efor
  elif strategy == 'all':
    for o in arr:
      h.update(str(o).encode())
    #efor
  #fi

  return h.hexdigest()
#edef

def uniq(lst, key=lambda x: x):
  """
  uniq: return unique elements in a list
  Input: lst : List of elements
         key : Function to perform on elements in lst (default is identity, example may be lambda x: str(x.attribute))
  Output: List of unique elements
  """
  U = {}
  for item in lst:
    k = key(item)
    if k not in U:
      U[k] = item
    #fi
  #efor
  return list(U.values())
#edef

def which(L, fn=lambda x: x, complement=False):
    """ Find the index of all elements in L that satisfy a certain criteria fn (identity by default). Complement returns all that satisfy opposite criteria. """
    if complement:
        wfn = lambda v: not(fn(v))
    else:
        wfn = fn
    #fi
    return [ idx for (idx, value) in enumerate(lst) if wfn(value) ]
#edef

def pairwise(lst, fun=lambda x, y: (x,y), symmetric=False, key=None):
    """
    Perform a function for each pair of elements in a list.
    lst : A list (or dictionary) of items
    fun : A function that takes two arguments A and B.
    symmetric : If the function is symmetric or not
    key: For each element in lst, perform this function on it first.
    
    output: ndarray of results of the function
    """
    keys = None
    if isinstance(lst, dict):
        keys = lst.keys()
        lst  = lst.values()
    #fi
    
    if key is not None:
        lst = [ [ key(e) for e in g] for g in lst ]
    #fi
    
    nelem = len(lst)
    
    pairwise = np.zeros((nelem, nelem))

    for set_i in range(nelem):
        for set_j in range(set_i, nelem):
            pairwise[set_i,set_j] = fun(lst[set_i], lst[set_j])
            pairwise[set_j,set_i] = pairwise[set_i,set_j] if symmetric else fun(lst[set_j], lst[set_i])
        #efor
    #efor
    return pairwise
#edef

def overlap(lst, key=None):
    """
    overlap: Return number of overlapping elements in each pair of sets
    Input: lst : A list of lists of elements
           key : Function to perform on each element in each list
    Output: A 2d numpy array of overlap counts
    """
    return pairwise(lst, fun=lambda A, B: len(set(A) & set(B)), symmetric=True, key=key).astype(int)
#edef

def jaccard(lst, key=None):
    """
    jaccard: Return jaccard index between each pair of sets
    Input: lst : A list of lists of elements
           key : Function to perform on each element in each list
    Output: A 2d numpy array of jaccard index
    """
    return pairwise(lst, fun=lambda A,B: float(len(set(A) & set(B))) / float(len(set(A) | set(B))), symmetric=True, key=key)
#edef


def chunks(l, m=None, n=None):
    """
    chunks: Yield successive m-sized chunks from l, or n (len(l)/n)-sized chunks from l
    inputs: l : List
            m : Maximum chunksize. stepsize -> m
            n : Number of chunks. stepsize -> math.ceil(len(l)/n)
      
      chunks the list l into stepsizes based on m or n. If neither m or n are specified, then stepsize will be 1
    Output: A generater of chunks
    """
        
    step = math.ceil(len(l)/n) if n is not None else m if m is not None else 1
        
    for i in range(0, len(l), step):
        yield l[i:i + step]
    #efor
#edef

def argrank(r, p=0):
    """
    argrank: Return the index of the first,second,nth value in a ranked list.
    
    parameters:
    -----------
    r: list[sortable value]
        A list of sortable values (numeric, strings, object)
    p: which position to take. Default is 0, meaning the highest value.
       1 would be the second highest value.
    """
    if isinstance(r, dict):
        return sorted(r.keys(), key=lambda x: r[x])[-1-p]
    #fi
    return sorted(enumerate(r), key=lambda x: x[1])[-1-p][0]
#edef
