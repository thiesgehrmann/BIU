"""
  Tools to analyze covariates.
   * expand_categorical
   * order_categories
   * dummy
   * cramers_corrected_stat 
   * associate_pair
   * associate
"""

from .. import utils
from .. import processing
from .. import stats

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")
sm = utils.py.loadExternalModule('statsmodels.api')
sstats = utils.py.loadExternalModule('scipy.stats')
sns    = utils.py.loadExternalModule('seaborn')

from collections import namedtuple
assoc_result = namedtuple('associationResult', ['method', 'statistic', 'pvalue'])

def expand_categorical(cov):
    """
    Expand categorical covariates into 
    Input:
      cov: Pandas DataFrame. covariates are annotated as categorical (with cov.col.astype('category') )
    Output:
      Pandas DataFrame, with categorical covariates removed, and new numerical covariates to replace it.
      e.g. original column: [ A,A,B,C ]:
      new columns:          [ 1,1,0,0 ]
                            [ 0,0,1,0 ]
                            [ 0,0,0,1 ]
    """
    E = cov.copy()
    for col in E.columns:
        if E[col].dtype.name != 'category':
            continue
        #fi
        values = sorted(E[col].drop_duplicates().values)
        for value in values:
            arr = np.zeros(len(E[col]))
            arr[np.where(E[col].values == value)] = 1
            E['%s_%s' % (col, str(value))] = arr
        #efor
        E = E.drop(columns=[col])
    #efor
    return E
#edef

def order_categories(categories, cont=None, statistic=np.mean):
    """
    Order categories based on a statistic of a continuous variable associated with the category
    Inputs:
        categories: a list of length n of categorical labels for n objects
        cont: a list of n continuous values for n objects
        statistic: function to calculate statistic
    Outputs:
        for all n objects, an integer with the position in the order
    """
    
    if isinstance(categories, pd.Series):
        categories = categories.values
    #fi
    if isinstance(cont, pd.Series):
        cont = cont.values
    elif cont is None:
        cont = np.ones(len(categories))
    #fi
    cats = set(categories)
    order = { cat : statistic(cont[np.where(categories == cat)]) for cat in cats }
    order = sorted(order.keys(), key=lambda x: order[x])
    order = { c : i for (i,c) in enumerate(order) }
    
    return np.array([ order[c] for c in categories])
#edef

def dummy(categorical_var):
    """
    Given a list of categorical variables, dummy code them
    Input:
      categorical_var : list of variables
    Output:
      matrix of dummy coded variables.
      e.g. Input:   [ A,A,B,C ]

                      A B C 
           Output: [[ 1,0,0 ],
                    [ 1,0,0 ],
                    [ 0,1,0 ],
                    [ 0,0,1 ]]
      
    """
    if not all([isinstance(cv, int) for cv in categorical_var ]):
        categorical_var = order_categories(categorical_var)
    d = np.zeros((len(categorical_var), max(categorical_var)+1))
    d[list(zip(*enumerate(categorical_var)))] = 1
    return d
#edef

    
def cramers_corrected_stat(confusion_matrix):
    """ calculate Cramers V statistic for categorial-categorial association.
        uses correction from Bergsma and Wicher, 
        Journal of the Korean Statistical Society 42 (2013): 323-328
        Taken from: https://stackoverflow.com/a/39266194
    """
    chi2 = sstats.chi2_contingency(confusion_matrix)[0]
    n = confusion_matrix.sum()
    phi2 = chi2/n
    r,k = confusion_matrix.shape
    phi2corr = max(0, phi2 - ((k-1)*(r-1))/(n-1))    
    rcorr = r - ((r-1)**2)/(n-1)
    kcorr = k - ((k-1)**2)/(n-1)
    statistic = np.sqrt(phi2corr / min( (kcorr-1), (rcorr-1)))
    pvalue = sstats.chi2.pdf(statistic, n)
    return assoc_result('cramers_v', statistic, pvalue)
#edef

def associate_pair(X, Y):
    """
    Determine if there is a significant association between two variables.
    For numeric variables it performs a regression.
    For categorical vs numeric variables it performs a regression (essentially an ANOVA)
    For categorical variables it performs a cramers V test
    
    It returns a p-value, and an effect size.
    """
    def cat_cat(d1, d2):
        d1 = order_categories(d1)
        d2 = order_categories(d2)
        cm = np.zeros((len(np.unique(d1)), len(np.unique(d2))))
        for i,j in zip(d1, d2):
            cm[int(i),int(j)] += 1
        #efor
        return cramers_corrected_stat(cm)
    #edef
    def cat_num(cat, num):
        #plt.scatter(order_categories(cat), num)
        cat = dummy(cat)
        cat = sm.add_constant(cat)
        num = num - np.mean(num) / np.std(num)
        ols_res = sm.OLS(num, cat).fit()
        effect, pvalue = min(zip(ols_res.params[1:], ols_res.pvalues[1:]), key=lambda x: x[1])
        return assoc_result('anova', effect, pvalue)
    #edef
        
    def num_cat(num, cat):
        return cat_num(cat, num)
    #edef
    
    def num_num(X, Y):
        X = (X - np.mean(X)) / np.std(X)
        Y = (Y - np.mean(Y)) / np.std(Y)
        X = sm.add_constant(X)
        ols_res = sm.OLS(Y,X).fit()
        return assoc_result('linreg', ols_res.params[1], ols_res.pvalues[1])
    #edef
    
    assocs = { (True, True) : cat_cat,
               (True, False) : cat_num,
               (False, True) : num_cat,
               (False, False) : num_num}
    X_no_na = X[~(pd.isna(X) | pd.isna(Y))]
    Y_no_na = Y[~(pd.isna(X) | pd.isna(Y))]
    
    return assocs[(X.dtype.name == 'category',Y.dtype.name == 'category')](X_no_na, Y_no_na)
#edef

def associate(covariates, data=None, nc=6, plot=False, ax=None, correctionType='fdr', **kwargs):
    """
    Correlate a matrix of covariates with itself, or with the first nc PCs of a data matrix
    Inputs:
      covariates: A dataframe of covariates (rows are data points, columns are covariates)
      data: A dataframe of data (rows are datapoints, columns are measurements)
      effect: Report the effect size instead of the pvalue
      nc : The number of components to inspect
      plot: Boolean. Plot a diagram of results if True
      ax: Matplotlib axis. Plot onto this axis, if None, generate own one.
      correctionType: The Multiple Testing Correction method to use (see biu.stats.correction.correct). You can specify None to not correct.
      **kwargs: Optional arguments for biu.stats.correction.correct
    Output:
        E : A dataframe of effects, and
        P : A dataframe of (corrected) p-values
    """
    E = None
    P = None
    if data is None:
        emat = np.zeros((len(covariates.columns), len(covariates.columns)))
        pmat = np.zeros((len(covariates.columns), len(covariates.columns)))
        for i, cov_i in enumerate(covariates.columns):
            for j,cov_j in enumerate(covariates):
                print('\r%s - %s %s' % (cov_i, cov_j, ''.join([' ']*30)), end='')
                res = associate_pair(covariates[cov_i].values, covariates[cov_j].values)
                emat[i,j] = res.statistic
                pmat[i,j] = res.pvalue
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=covariates.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=covariates.columns)

    else:
        pc_data = processing.dataframe.pca(data, nc=nc)
        emat = np.zeros((len(covariates.columns), nc))
        pmat = np.zeros((len(covariates.columns), nc))
        for i, cov_i in enumerate(covariates.columns):
            for j, comp_j in enumerate(pc_data.columns):
                print('\r%s - %s %s' % (cov_i, comp_j, ''.join([' ']*30)), end='')
                res = associate_pair(pc_data[comp_j].values, covariates[cov_i].values)
                emat[i,j] = res.statistic
                pmat[i,j] = res.pvalue
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=pc_data.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=pc_data.columns)
    #fi

    if correctionType is not None:
        P = stats.correction.correct(P, correctionType=correctionType, **kwargs)
    #fi

    if plot:
        if ax is None:
            fig, axes = utils.figure.subplots(nrows=1, ncols=1, figsize=np.array(P.shape[::-1])/4)
            ax = axes[0]
        #fi
        sns.heatmap(E, ax=ax, center=0, cmap="PiYG", fmt='',
                   annot=P.applymap(lambda x: '*' if x < 0.05 else ''))
    #fi

    return E, P
#edef
