from .. import utils
from . import matrix

skd = utils.py.loadExternalModule('sklearn.decomposition')
skm = utils.py.loadExternalModule('sklearn.manifold')
pd  = utils.py.loadExternalModule('pandas')

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

def reorder(df, **kwargs):
  """
  Reorder the columns/rows of a dataframe
  See biu.ops.matrix.reorder for information on this function
  """
  return matrix.reorder(df, **kwargs)
#edef
