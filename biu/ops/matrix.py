
from .. import utils

fastcluster = utils.py.loadExternalModule('fastcluster')
np = utils.py.loadExternalModule('numpy')
pd = utils.py.loadExternalModule('pandas')
ssdist = utils.py.loadExternalModule('scipy.spatial.distance')
sstats = utils.py.loadExternalModule('scipy.stats')

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

####################################################################

def corrcoef(matrix, axis=0, method='pearson'):
    """
    Calculate pearsons correlation coefficient for a numpy matrix.
    Returns a matrix of correlation coefficients, and estimated p-values per correlation
    
    parameters:
    -----------
    matrix: A numpy matrix to calculate correlations for
    axis: Which axis to perform the operation on.
          0 : rows
          1 : columns
    method=['pearson', 'spearman' ] or callable which must return (r-correlation, and p-value) tuple

    
    returns:
    --------
    r: Correlation coefficients
    p: p-values
    
    Modified from: https://stackoverflow.com/a/24547964
    """
    
    def _corrcoef_pearson(matrix, axis):
        from scipy.special import betainc

        if axis == 1:
            matrix = matrix.transpose()
        #fi

        r = np.corrcoef(matrix)
        rf = r[np.triu_indices(r.shape[0], 1)]
        df = matrix.shape[1] - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = betainc(0.5 * df, 0.5, df / (df + ts))
        p = np.zeros(shape=r.shape)
        p[np.triu_indices(p.shape[0], 1)] = pf
        p[np.tril_indices(p.shape[0], -1)] = pf
        p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])

        return r, p
    #edef
    
    def _corrcoef_loop(matrix, axis, method):
        
        if axis == 1:
            matrix = matrix.transpose()
        #fi
        
        rows, cols = matrix.shape[0], matrix.shape[1]
        r = np.ones(shape=(rows, rows))
        p = np.ones(shape=(rows, rows))
        for i in range(rows):
            for j in range(i+1, rows):
                r_, p_ = method(matrix[i], matrix[j])
                r[i, j] = r[j, i] = r_
                p[i, j] = p[j, i] = p_
            #efor
        #efor
        return r, p
    #edef
    
    method = method.lower()
    
    if method == 'pearson':
        return _corrcoef_pearson(matrix, axis)
    elif method == 'spearman':
        res = sstats.spearmanr(matrix, axis=1-axis)
        return res.correlation, res.pvalue
    elif hasattr(method, '__call__'):
        return _corrcoef_loop(matrix, axis, method)
    else:
        raise ValueError("Don't know what to do with method '%s'" % method)
    #fi
#edef

####################################################################

def corrcoef_between(matrix1, matrix2, axis=0, method='pearson'):
    """
    Calculate the correlation between columns/rows in two different matrices.
    The opposite axis must have the same dimension in the two matrices.
    
    parameters:
    -----------
    matrix1: A numpy matrix 
    matrix2: A numpy matrix
    axis: Which axis to perform the operation on.
          0 : rows
          1 : columns
    method=['pearson', 'spearman' ] or callable which must return (r-correlation, and p-value) tuple
    
    returns:
    --------
    r: Correlation coefficients
    p: p-values
    
    Modified from: https://stackoverflow.com/a/24547964
    """
    
    # NOTE: I might want to think about concatenating the matrices for the pearson/spearman methods
    #       and then extracting the relevant columns/rows
    
    def _corrcoef_between_callable(matrix1, matrix2, axis, method):
        rows = matrix1.shape[axis]
        cols = matrix2.shape[axis]

        if matrix1.shape[1-axis] != matrix2.shape[1-axis]:
            raise ValueError("Cannot calculate correlation along this axis.")
        #fi

        r = np.ones(shape=(rows, cols))
        p = np.ones(shape=(rows, cols))
        for i in range(rows):
            dat1 = matrix1[i,:] if (axis == 0) else matrix1[:,i]
            for j in range(cols):
                dat2 = matrix2[j,:] if (axis == 0) else matrix2[:,j]
                r_, p_ = method(dat1, dat2)
                r[i, j] = r_
                p[i, j] = p_
            #efor
        #efor
        return r, p
    #edef
    
    def _corrcoef_spearman(matrix1, matrix2, axis):
        r, p = sstats.spearmanr(matrix1, matrix2, 1-axis)
        
        if axis == 1:
            r = r[:matrix1.shape[axis], :][:, matrix1.shape[axis]:]
            p = p[:matrix1.shape[axis], :][:, matrix1.shape[axis]:]
        else:
            r = r[:, :matrix1.shape[axis]][matrix1.shape[axis]:, :]
            p = p[:, :matrix1.shape[axis]][matrix1.shape[axis]:, :]
        #fi

        return r, p
    #edef
        
        

    method = method.lower()
    
    if method == 'pearson':
        return _corrcoef_between_callable(matrix1, matrix2, axis=axis, method=sstats.pearsonr)
    elif method == 'spearman':
        return _corrcoef_spearman(matrix1, matrix2, axis)
    elif method == 'spearman_slow':
        return _corrcoef_between_callable(matrix1, matrix2, axis=axis, method=sstats.spearmanr)
    elif hasattr(method, '__call__'):
        return _corrcoef_between_callable(matrix1, matrix2, axis=axis, method=method)
    else:
        raise ValueError("Don't know what to do with method '%s'" % method)
    #fi

#edef
