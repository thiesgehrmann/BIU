from .. import utils
from . import matrix

from . import series

skd = utils.py.loadExternalModule('sklearn.decomposition')
skm = utils.py.loadExternalModule('sklearn.manifold')
pd  = utils.py.loadExternalModule('pandas')
np  = utils.py.loadExternalModule('numpy')

#########################################################################

def detect_categorical(df, ret_types=False):
    """
    Detects and sets categorical columns in a dataframe.
    Inputs:
        df: A pandas DataFrame
        ret_types: Return a dictionary of column: [category|numeric]
    Outputs:
        A copy of the dataframe, in which columns detected to be categorical have the appropriate dtype set ('category')
        Numeric columns (columns with numbers) are coerced to be numeric floats
        
        if ret_types is True:
         (dataframe, dict)
         
    Note: This function detects if a column is numeric by testing if each non-NA cell contains an
          integer/float, or an integer-like string or a float-like string. Otherwise, it is
          set as categorical.
    """
    new_df = df.copy()
    types = {}
    for col in df.columns:
        S = df[col]
        if series.is_categorical(S):
            new_df[col] = series.cast_category(S)
            types[col] = 'category'
        else:
            new_df[col] = series.cast_float(S)
            types[col] = 'numeric'
        #fi
    #efor
    
    if ret_types:
        return new_df, types
    else:
        return new_df
    #fi
#efor

#########################################################################

def pca(df, nc=2, ret_fit=False, **kwargs):
    """
    Perform PCA on a dataframe
    Inputs:
        df: DataFrame
        nc: The number of components
        ret_fit: Boolean. Return the fit, together with the dataframe.
        **kwargs: Options for Scikit-learn PCA decomposition
    Outputs:
        if ret_fit is False:
            DataFrame with index preserved, columns are PCA_1, PCA_... PCA_n
        else:
            Same dataframe as above, plus the fit object
    """
    fit = skd.PCA(nc, **kwargs).fit(df.values)
    DF  = pd.DataFrame(fit.transform(df.values), index=df.index, columns=[ 'PCA_%d' % (i+1) for i in range(nc)])
    if ret_fit:
        return DF, fit
    else:
        return DF
    #fi
#edef

#########################################################################

def tsne(data, **kwargs):
    """
    Perform TSNE on a matrix/dataframe
    Inputs:
        data: Dataframe/matrix/ndarray
        **kwargs: Options for TSNE
    Outputs:
        DataFrame/matrix of TSNE reduction
    """
    trans = skm.TSNE(**kwargs).fit_transform(data)
    trans = pd.DataFrame(trans, index=data.index, columns=[ 'TSNE_%d' % (i+1) for i in range(trans.shape[1])])
    return trans

#edef

#########################################################################

def reorder(df, **kwargs):
  """
  Reorder the columns/rows of a dataframe
  See biu.ops.matrix.reorder for information on this function
  """
  return matrix.reorder(df, **kwargs)
#edef

#########################################################################

def multimax(A, B, *pargs):
    """
    Determine the maximum value per cell in a group of dimension-matched dataframes
    
    Parameters:
    -----------
    A, B: pandas.DataFrame
        
    *pargs: More DataFrames
    
    Output:
    -------
    pandas.DataFrame, with index and columns of input A
    """
    I = [ A, B ] + list(pargs)
    M = matrix.multimax(*[ a.values for a in I])
    return pd.DataFrame(M, columns=A.columns, index=A.index)
#edef

#########################################################################

def flat(df, fields=None, sep=None, ignore_unequal_rows=False):
    """
    Flatten a Pandas DataFrame, given columns with lists in them
    Parameters:
    -----------
    
    df:
        Pandas DataFrame
    fields:
        The fields to flatten
    sep:
        If the fields are a list in string format, you can specify the delimiter here
    ignore_uneqial_rows:
        If the lists are not exactly the same length, you can choose to ignore the error induced by this here.
        
    Output:
        Pandas DataFrame with flattened columns
    
    example 1:
    ---------------------
    
    x = idx | f1 | f2    | f3
          0 |  1 | [0,1] | [2,3]
          1 | 34 | [8,9] | [7,6]
          
    flat(x, ['f1','f2']) =
        idx | f1 | f2 | f3
          0 |  1 |  0 |  2
          0 |  1 |  1 |  3
          1 | 34 |  8 |  7
          1 | 34 |  9 |  6
          
    example 2:
    ----------------------
    
    x = idx | f1    | f2
          0 | [0,1] | [2,3]
          1 | [8,9] | [7,6]
          
    flat(x) =
        idx | f1 | f2
          0 |  0 |  2
          0 |  1 |  3
          1 |  8 |  7
          1 |  9 |  6
    """
    
    if fields is None:
        fields = list(df.columns)
    #fi
    
    non_grouped_fields = [ f for f in df.columns if f not in fields ]

    grouped = df.copy()
    
    if sep is not None:
        for f in fields:
            grouped[f] = grouped[f].apply(lambda x: x.split(sep))
        #efor
    #fi
    
    for f in fields:
        if not hasattr(grouped[f].iloc[0], '__iter__'):
            raise ValueError("Field %s doesn't is not iterable")
        #fi
    #efor
    
    lengths_match = grouped[fields].applymap(len).apply(lambda x: np.all(x == x[0]), axis=1)
    
    if (not lengths_match.all()) and (not ignore_unequal_rows):
        raise ValueError("Some rows don't have the same length in all columns. Set ignore_unequal_rows=True if you want to ignore.")
    #fi

    flattened = []
    indexes   = []
    for (i, row) in grouped.iterrows():
        for values in zip(*row[fields]):
            flattened.append(list(values) + list(row[non_grouped_fields].values))
            indexes.append(i)
        #efor
    #efor

    flattened = pd.DataFrame(flattened, index=indexes, columns=fields + non_grouped_fields)
    
    return flattened
#edef

#########################################################################

def corrcoef(dataframe, axis=0, method='pearson'):
    """
    Calculate pearsons correlation coefficient for the columns in a dataframe.
    Returns a matrix of correlation coefficients, and estimated p-values per correlation
    
    parameters:
    -----------
    dataframe: A pandas dataframe to calculate correlations for
    axis: Which axis to perform the operation on.
          0 : rows
          1 : columns
    method=['pearson', 'spearman' ] or callable which must return (r-correlation, and p-value) tuple
    
    returns:
    --------
    r: Correlation coefficients
    p: p-values
    
    Found at: https://stackoverflow.com/a/24547964
    """
    
    r, p = matrix.corrcoef(dataframe.values, axis, method=method)
    idx = dataframe.index if axis==0 else dataframe.columns
    r = pd.DataFrame(r, columns=idx, index=idx)
    p = pd.DataFrame(p, columns=idx, index=idx)
    
    return r, p
#edef

####################################################################

def corrcoef_between(dataframe1, dataframe2, axis=0, method='pearson'):
    """
    Calculate the correlation coefficient between columns/rows in two different matrices.
    The opposite axis must have the same dimension in the two matrices.
    
    parameters:
    -----------
    dataframe1: A pandas dataframe 
    dataframe2: A pandas dataframe
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
    
    r, p = matrix.corrcoef_between(dataframe1.values, dataframe2.values, axis=axis, method=method)
    idx = dataframe1.index if (axis == 0) else dataframe1.columns
    col = dataframe2.index if (axis == 0) else dataframe2.columns

    r = pd.DataFrame(r, columns=col, index=idx)
    p = pd.DataFrame(p, columns=col, index=idx)
    
    return r, p
#edef