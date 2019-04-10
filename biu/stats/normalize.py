from .. import utils


np = utils.py.loadExternalModule('numpy')
pd = utils.py.loadExternalModule('pandas')
sstats = utils.py.loadExternalModule('scipy', 'stats')

from .common import *

#################################################################################

def rin(D, axis=0):
    """
    Perform a Rank Inverse Normal normalization on a dataframe, skipping the NA values.
    
    parameters:
    D: Pandas DataFrame.
    axis : {0 or 'index', 1 or 'columns'}, default 0
        Axis along which the function is applied:

        * 0 or 'index': apply function to each column.
        * 1 or 'columns': apply function to each row.
        
    returns:
    Pandas DataFrame
    """
    def single_rin(x):
        y = x.copy()
        nas = pd.isna(x)
        r = sstats.rankdata(x[~nas])
        r = (r-0.5) / len(x)
        r = q_(r)
        y[~nas] = r
        return y
    #edef
    
    # Perfor the action on the columns, but for some reason I need to do a transpose first...
    if (axis == 1) or (axis == 'columns'):
        return D.transpose().apply(single_rin, axis=0).transpose()
    else:
        return D.apply(single_rin, axis=0)
    #fi
#edef

#################################################################################

def zscore(df, axis=0):
    """
    Perform a zscore normalization on a dataframe, skipping the NA values.
    
    parameters:
    df: Pandas DataFrame.
    axis : {0 or 'index', 1 or 'columns'}, default 0
        Axis along which the function is applied:

        * 0 or 'index': apply function to each column.
        * 1 or 'columns': apply function to each row.
        
    returns:
    Pandas DataFrame
    """
    def single_zscore(A):
        B = A.copy()
        nas = pd.isna(A)
        B[~nas] = sstats.zscore(A[~nas])
        return B
    #edef
    
    return df.apply(single_zscore, axis)
#edef

#################################################################################