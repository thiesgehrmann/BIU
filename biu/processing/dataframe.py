from .. import utils
from . import matrix

skd = utils.py.loadExternalModule('sklearn.decomposition')
skm = utils.py.loadExternalModule('sklearn.manifold')
pd = utils.py.loadExternalModule('pandas')

def pca(df, nc=2, **kwargs):
    """
    Perform PCA on a dataframe
    Inputs:
        df: DataFrame
        nc: The number of components
        **kwargs: Options for Scikit-learn PCA decomposition
    Outputs:
        DataFrame with index preserved, columns are PCA_1, PCA_... PCA_n
    """
    return pd.DataFrame(skd.PCA(nc, **kwargs).fit(df.values).transform(df.values), index=df.index, columns=[ 'PCA_%d' % (i+1) for i in range(nc)])
#edef

from sklearn.manifold import TSNE

def tsne(data, **kwargs):
    """
    Perform TSNE on a matrix/dataframe
    Inputs:
        data: Dataframe/matrix/ndarray
        **kwargs: Options for TSNE
    Outputs:
        DataFrame/matrix of TSNE reduction
    """
    t = skm.TSNE().fit_transform(bloodExpr[diffExGenes])
    if isinstance(data, pd.DataFrame):
        t = pd.DataFrame(t, index=data.index, columns=[ 'TSNE_%d' % (i+1) for i in range(x.shape[1])])
    #fi
    return t
#edef
        

def reorder(df, **kwargs):
  """
  Reorder the columns/rows of a dataframe
  See biu.processing.matrix.reorder for information on this function
  """
  return matrix.reorder(df, **kwargs)
#edef
