from .. import utils
from .. import math

from math import ceil

np     = utils.py.loadExternalModule("numpy")
linalg = utils.py.loadExternalModule("scipy.linalg")
sm     = utils.py.loadExternalModule('statsmodels.api') # For some reason I need to import this before I can import the line below...
smf    = utils.py.loadExternalModule('statsmodels.formula', 'api')

###############################################################################

def lowess(x, y, f=2. / 3., iter=3, function=False):
    """lowess(x, y, f=2./3., iter=3) -> yest

    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.

    For more information, see
    
    William S. Cleveland: "Robust locally weighted regression and smoothing
    scatterplots", Journal of the American Statistical Association, December 1979,
    volume 74, number 368, pp. 829-836.
    
    William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
    approach to regression analysis by local fitting", Journal of the American
    Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
    
    # Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
    #
    # License: BSD (3-clause)

    function : Return a piece wise linearly interpolated function of the lowess curve.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]
        #efor

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2
    #efor

    if function:
        x, yest = zip(*sorted(zip(x, yest), key=lambda x: x[0])) 
        return lambda position: math.interpolation.linearInterpolation(x, yest, position, isSorted=True)
    #fi

    return yest
#edef

###############################################################################

def linear_residuals(formula, Y, covariates, intercept=True):
    """
    Regress out a set of covariates from a set of measurements.
    For each column in Y_columns, regress out the variables in covariates, given the formula.
    (These are done independently per measurement...)

    Inputs:
        formula:    R style formula string (i.e. ~ 0 + cov1 + cov2)
                    If None, then the linear summation of all values in covariates is used.
                    C() around categorical parameters are not necessary, if they are encoded as categorical in the covariates dataframe.
        Y:          Pandas Data Frame of measurements
                    Columns of values (e.g. expression.) Columns are individual measures, and rows are samples
                    e.g. rows are individuals and columns are gene expression
        covariates: Pandas DataFrame of covariates
                    Columns of value (e.g. BMI) Columns are covariates, and rows are sample
                    The index must correspond to that in Y_columns
        intercept:  If no formula is specified, then add an intercept/bias/mean term to the generated formula.
    Output:
        Y_resids
    """
    
    column_rename = { c: 'meas_%d' % (i+1) for i,c in enumerate(Y.columns) }
    column_dename = {v: k for k, v in column_rename.items()}

    resids = Y.copy().rename(columns=column_rename)
    data   = covariates.join(resids)

    if formula is None:
        formula = "~ %s%s" % ( ('' if intercept else '0 + '), ' + '.join(covariates.columns))
    elif '~' not in formula:
        formula = '~ %s' % formula
    #fi

    formula = formula.split('~')[1].strip()

    for cov in covariates.columns:
        if (covariates[cov].dtype.name == 'category') and ('C(%s)' % cov not in formula) and (cov in formula):
            formula = formula.replace(cov, 'C(%s)' % cov)
        #fi
    #efor

    for i, gene in enumerate(resids.columns):
        gene_formula = "%s ~ %s" % (gene, formula)
        resids[gene] = smf.ols(formula=gene_formula, data=data).fit().resid
        print('\r(%d) [ %s ~ %s ]%s' % (i+1, column_dename[gene], formula, ' '*50), end="")
    #efor
    
    resids = resids.rename(columns=column_dename)

    return resids
#edef
