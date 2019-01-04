
from .. import utils

fastcluster = utils.py.loadExternalModule('fastcluster')
np = utils.py.loadExternalModule('numpy')
pd = utils.py.loadExternalModule('pandas')
ssdist = utils.py.loadExternalModule('scipy.spatial.distance')

def order(M, distance="correlation", method="single"):
  '''
   input:
     - M is a matrix
     - distance is the desired distance metric used by pdist
     - method is the linkage method
   output:
     - order implied by the dendrogram produced
  '''

  def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
      return [cur_index]
    else:
      left = int(Z[cur_index-N,0])
      right = int(Z[cur_index-N,1])
      return (seriation(Z,N,left) + seriation(Z,N,right))
    #fi
  #edef
    
  def computeSerialMatrix(flat_dist_mat, N, method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    res_linkage = fastcluster.linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    return res_order
  #edef

  M = np.matrix(M)
  D = ssdist.pdist(M, distance)
  D[np.isnan(D)] = -1
  return computeSerialMatrix(D, len(M), method=method)
#edef

####################################################################

def reorder(matrix, symmetric=False, by=None, **kwargs):
    """
    Reorder an ndarray, matrix or a dataframe
    Inputs:
        matrix: pandas DataFrame, numpy ndArray or numpy matrix
        symmetric: If it is a symmetric matrix, reorder the columns as well as the rows
        by : Reorder this matrix by the ordering of this matrix, instead
        **kwargs: The options for biu.ops.matrix.order (distance, method etc.)
    Output:
        Depending on input:
            * pandas DataFrame
            * numpy ndArray
            * numpy matrix
    """
    if isinstance(matrix,pd.DataFrame):
        ordered = order(matrix if by is None else by, **kwargs)
        if symmetric:
            return matrix.iloc[ordered][[matrix.columns[i] for i in ordered]]
        else:
            return matrix.iloc[ordered]
        #fi
    elif isinstance(matrix, np.ndarray) or isinstance(matrix, np.matrix):
        ordered = order(matrix if by is None else by, **kwargs)
        if symmetric:
            return matrix[order,ordered]
        else:
            return matrix[order,:]
        #fi
    else:
        raise TypeError
    #fi
#edef