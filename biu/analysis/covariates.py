"""
  Tools to analyze covariates.
   * detect_categorical
   * expand_categorical
   * order_categories
   * dummy
   * cramers_corrected_stat 
   * associate_pair
   * associate

"""

from .. import utils
from .. import ops
from .. import stats

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")
sm = utils.py.loadExternalModule('statsmodels.api')
sstats = utils.py.loadExternalModule('scipy.stats')
sns    = utils.py.loadExternalModule('seaborn')

from collections import namedtuple
assoc_result = namedtuple('associationResult', ['method', 'statistic', 'pvalue'])

##############################################################################################

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
    return ops.dataframe.detect_categorical(df, ret_types)
#edef

##############################################################################################

def expand_categorical(cov, sep='_', as_bool=False, exclude=[]):
    """
    Expand categorical covariates into 
    Input:
      cov: Pandas DataFrame. covariates are annotated as categorical (with cov.col.astype('category') )
      sep: String. In the new column names, separate the column name and the value with this charachter.
              e.g. column name "happy" with values [ yes, no, maybe] and sep '_' will produce
              happy_yes, happy_no, happy_maybe
      as_bool: Make dummy variables with True/False instead of 1/0
      exclude: List. Exclude these columns from the transformation.
    Output:
      Pandas DataFrame, with categorical covariates removed, and new numerical covariates to replace it.
      e.g. original column: c [ A,A,B,C ]:
      new columns:          c_A [ 1,1,0,0 ]
                            c_B [ 0,0,1,0 ]
                            c_C [ 0,0,0,1 ]
    """
    E = cov.copy()
    for col in E.columns:
        if E[col].dtype.name != 'category':
            continue
        elif col in exclude:
            continue
        elif all([ isinstance(x, bool) for x in E[col] if not pd.isna(x) ]):
            continue
        #fi
        values = sorted(E[col].cat.remove_unused_categories().cat.categories.values)
        for value in values:
            arr = np.zeros(len(E[col]))
            arr[E[col].isna()] = None
            arr[np.where(E[col].values == value)] = 1
            arr = [ True if v == 1 else False if v == 0 else None for v in arr ] if as_bool else arr
            E['%s%s%s' % (col, sep, str(value))] = arr
            E['%s%s%s' % (col, sep, str(value))] = E['%s%s%s' % (col, sep, str(value))].astype('category')
        #efor
        E = E.drop(columns=[col])
    #efor
    return E
#edef

##############################################################################################

def order_categories(categories, cont=None, statistic=None):
    """
    Order categories based on a statistic of a continuous variable associated with the category
    Inputs:
        categories: a list of length n of categorical labels for n objects
        cont: a list of n continuous values for n objects
        statistic: function to calculate statistic (default is np.mean)
    Outputs:
        for all n objects, an integer with the position in the order
    """
    
    statistic = np.mean if statistic is None else statistic
    
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

##############################################################################################

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
    cv_orig = categorical_var
    if not all([isinstance(cv, int) for cv in categorical_var ]):
        categorical_var = order_categories(categorical_var)
    d = np.zeros((len(categorical_var), max(categorical_var)+1))
    d[tuple(zip(*enumerate(categorical_var)))] = True
    d[pd.isna(cv_orig),:] = None
    return d[:,sorted(set(categorical_var))]
#edef

##############################################################################################
    
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

##############################################################################################

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
        
        X = X[np.isfinite(X) & np.isfinite(Y)]
        Y = Y[np.isfinite(X) & np.isfinite(Y)]
        
        X = (X - np.mean(X)) / np.std(X)
        Y = (Y - np.mean(Y)) / np.std(Y)
        X = sm.add_constant(X)

        try: 
            ols_res = sm.OLS(Y,X).fit()
            return assoc_result('linreg', ols_res.params[1], ols_res.pvalues[1])
        except Exception as e:
            utils.msg.dbm("Encountered singular value. Cannot associate for this one. Returning no association.")
            return assoc_result('linreg_singular', 0, 1)
        #etry
    #edef
    
    assocs = { (True, True) : cat_cat,
               (True, False) : cat_num,
               (False, True) : num_cat,
               (False, False) : num_num}
    
    # Remove NA and INfs
    X_no_na = X[~pd.isna(X) & ~pd.isna(Y)]
    Y_no_na = Y[~pd.isna(X) & ~pd.isna(Y)]
    
    xtype = X.dtype.name == 'category'
    ytype = Y.dtype.name == 'category'
    
    if all([ isinstance(x, bool) for x in X_no_na ]):
        xtype = False
        X_no_na = X_no_na.astype(bool)
    #fi
    
    if all([ isinstance(y, bool) for y in Y_no_na ]):
        ytype = False
        Y_no_na = Y_no_na.astype(bool)
    #fi
                  
    return assocs[(xtype,ytype)](X_no_na, Y_no_na)

#edef

##############################################################################################

def associate(covariates, data=None, pca=True, nc=6, plot=False, ax=None, method='fdr', **kwargs):
    """
    Correlate a matrix of covariates with itself, or with the first nc PCs of a data matrix
    Inputs:
      covariates: A dataframe of covariates (rows are data points, columns are covariates)
      data: A dataframe of data (rows are datapoints, columns are measurements)
            If data is None, then each covariate will be associated with the other covariates
      effect: Report the effect size instead of the pvalue
      pca: Boolean. If true, perform PCA on data first. If False, do not.
      nc : The number of components to inspect
      plot: Boolean. Plot a diagram of results if True
      ax: Matplotlib axis. Plot onto this axis, if None, generate own one.
      method: The Multiple Testing Correction method to use (see biu.stats.p_adjust). You can specify None to not correct.
      **kwargs: Optional arguments for biu.stats.p_adjust
    Output:
        E : A dataframe of effects, and
        P : A dataframe of (corrected) p-values
    """
    E = None
    P = None
    
    covariates = expand_categorical(covariates, as_bool=True)
    if not pca:
        data = expand_categorical(data, as_bool=True)
    #fi
    
    if data is None:
        emat = np.zeros((len(covariates.columns), len(covariates.columns)))
        pmat = np.zeros((len(covariates.columns), len(covariates.columns)))
        for i, cov_i in enumerate(covariates.columns):
            for j,cov_j in enumerate(covariates):
                print('\r%s - %s %s' % (cov_i.replace('\n',' '), cov_j.replace('\n',' '), ''.join([' ']*30)), end='')
                try:
                    res = associate_pair(covariates[cov_i].values, covariates[cov_j].values)
                    emat[i,j] = res.statistic
                    pmat[i,j] = res.pvalue
                except ValueError:
                    emat[i,j] = 0
                    pmat[i,j] = 1
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=covariates.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=covariates.columns)

    else:
        data_fmt = data
        if pca:
            data_fmt = ops.dataframe.pca(data, nc=nc)
        #fi
        emat = np.zeros((len(covariates.columns), len(data_fmt.columns)))
        pmat = np.zeros((len(covariates.columns), len(data_fmt.columns)))
        for i, cov_i in enumerate(covariates.columns):
            for j, comp_j in enumerate(data_fmt.columns):
                print('\r%s - %s %s' % (cov_i, comp_j, ''.join([' ']*30)), end='')
                res = associate_pair(data_fmt[comp_j].values, covariates[cov_i].values)
                emat[i,j] = res.statistic
                pmat[i,j] = res.pvalue
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=data_fmt.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=data_fmt.columns)
    #fi

    if method is not None:
        P = stats.p_adjust(P, method=method, **kwargs)
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

##############################################################################################


def _correlation(X, Y):
    valid = ~(pd.isna(X) | pd.isna(Y))
    x = X[valid].astype(float).values
    y = Y[valid].astype(float).values
    r = np.corrcoef(x, y)[0,1]
    p = None
    if r == 1:
        p = 0
    else:
        n = sum(valid)
        t = np.abs(r) * np.sqrt(n-2) / np.sqrt(1-r**2)
        p = sstats.t.cdf(t, n-2)
    #fi
    return assoc_result('pearson', r, p)
#edef

def correlate(covariates, data=None, pca=True, nc=6, plot=False, ax=None, method='fdr', **kwargs):
    """
    Correlate a matrix of covariates with itself, or with the first nc PCs of a data matrix
    Inputs:
      covariates: A dataframe of covariates (rows are data points, columns are covariates)
      data: A dataframe of data (rows are datapoints, columns are measurements)
            If data is None, then each covariate will be associated with the other covariates
      effect: Report the effect size instead of the pvalue
      pca: Boolean. If true, perform PCA on data first. If False, do not.
      nc : The number of components to inspect
      plot: Boolean. Plot a diagram of results if True
      ax: Matplotlib axis. Plot onto this axis, if None, generate own one.
      method: The Multiple Testing Correction method to use (see biu.stats.p_adjust). You can specify None to not correct.
      **kwargs: Optional arguments for biu.stats.p_adjust
    Output:
        E : A dataframe of effects, and
        P : A dataframe of (corrected) p-values
    """
    E = None
    P = None
    
    covariates = expand_categorical(covariates)
    if not pca:
        data = expand_categorical(data)
    #fi
    
    if data is None:
        emat = np.zeros((len(covariates.columns), len(covariates.columns)))
        pmat = np.zeros((len(covariates.columns), len(covariates.columns)))
        for i, cov_i in enumerate(covariates.columns):
            for j,cov_j in enumerate(covariates):
                print('\r%s - %s %s' % (cov_i, cov_j, ''.join([' ']*30)), end='')
                x = covariates[cov_j]
                y = covariates[cov_i]
                res = _correlation(x, y)
                emat[i,j] = res.statistic
                pmat[i,j] = res.pvalue
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=covariates.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=covariates.columns)

    else:
        data_fmt = data
        if pca:
            data_fmt = ops.dataframe.pca(data, nc=nc)
        #fi
        emat = np.zeros((len(covariates.columns), len(data_fmt.columns)))
        pmat = np.zeros((len(covariates.columns), len(data_fmt.columns)))
        for i, cov_i in enumerate(covariates.columns):
            for j, comp_j in enumerate(data_fmt.columns):
                print('\r%s - %s %s' % (cov_i, comp_j, ''.join([' ']*30)), end='')
                x = data_fmt[comp_j]
                y = covariates[cov_i]
                res = _correlation(x, y)
                
                emat[i,j] = res.statistic
                pmat[i,j] = res.pvalue
            #efor
        #efor
        E = pd.DataFrame(emat, index=covariates.columns, columns=data_fmt.columns)
        P = pd.DataFrame(pmat, index=covariates.columns, columns=data_fmt.columns)
    #fi

    if method is not None:
        P = stats.p_adjust(P, method=method, **kwargs)
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

##############################################################################################
