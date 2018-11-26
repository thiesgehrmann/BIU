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

from sklearn.manifold import TSNE

def tsne(data, ret_fit=False, **kwargs):
    """
    Perform TSNE on a matrix/dataframe
    Inputs:
        data: Dataframe/matrix/ndarray
        ret_fit: Boolean. Return the fit, together with the dataframe
        **kwargs: Options for TSNE
    Outputs:
        if ret_fit is False:
            DataFrame/matrix of TSNE reduction
        else:
            Same dataframe as above, plus the fit object
    """
    rel_df = data[diffExGenes]
    fit    = skm.TSNE(**kwargs).fit(rel_df)
    trans  = fit.transform(rel_df)
    if isinstance(data, pd.DataFrame):
        trans = pd.DataFrame(trans, index=data.index, columns=[ 'TSNE_%d' % (i+1) for i in range(trans.shape[1])])
    #fi
    if ret_fit:
        return trans, fit
    else:
        return trans
    #fi
#edef
        

def reorder(df, **kwargs):
  """
  Reorder the columns/rows of a dataframe
  See biu.ops.matrix.reorder for information on this function
  """
  return matrix.reorder(df, **kwargs)
#edef
