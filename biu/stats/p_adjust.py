from .. import utils
from ..ops.array import *

np = utils.py.loadExternalModule("numpy")
pd = utils.py.loadExternalModule("pandas")

def p_adjust(p, method='fwer', n=None):
    """
    Adjust p values.
    This is a re-write of the R function p.adjust
    
    Inputs:
     - p: The array of pvalues
     - method: bonferroni, fwer
              holm
              hommel (not yet implemented
              hochberg
              benjamini_hochberg, fdr, fdr_bh
              banjamini_yekutieli, fdr_by
     - n: number of tests to correct for

    Outputs:
     - Corrected q-values
      
    """
    porig = p

    
    if isinstance(p, pd.DataFrame):
        p = p.values.reshape(np.product(p.shape))
    elif isinstance(p, np.ndarray) or isinstance(p, np.matrix):
        p = p.reshape(np.product(p.shape))
    else:
        p = np.array(p)
    #fi
    
    if n is None:
        n = len(p)
    #fi
    
    methods = { 'bonferroni': bonferroni, 'fwer' : bonferroni,
                'holm' : holm,
                'hommel' : hommel,
                'hochberg' : hochberg,
                'benjamini_hochberg' : BH, "bh" : BH, 'fdr' : BH, 'fdr_bh' : BH,
                'bejamini_yekutieli' : BY, 'by' : BY, 'fdr_by' : BY}
    
    method = method.lower()
    if method not in methods:
        method = 'fwer'
    #fi
    
    lp = len(p)
    if n < lp:
        n = lp
    #fi
    if (n == 2) & (method == 'hommel'):
        method = 'hochberg'
    #fi
    
    q = np.ones(lp)
    q[~np.isnan(p)] = methods[method](p[~np.isnan(p)], n)
    
    return q
#edef
    
    

def bonferroni(p, n):
    """
    Bonferroni correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
        bonferroni = pmin(1, n * p),"""
    return pmin(1, n * p)
#edef

def holm(p, n):
    """
    Holm correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
        holm = {
        i <- seq_len(lp)
        o <- order(p)
        ro <- order(o)
        pmin(1, cummax( (n - i + 1L) * p[o] ))[ro]
        },
    """
    lp = len(p)
    i = np.array(list(range(lp)))+1
    o = order(p)
    ro = order(o)
    return pmin(1, cummax((n - I + 1) * p[o]))[ro]
#edef

def hommel(p, n):
    """
    Hommel correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
        hommel = { ## needs n-1 >= 2 in for() below
        if(n > lp) p <- c(p, rep.int(1, n-lp))
        i <- seq_len(n)
        o <- order(p)
        p <- p[o]
        ro <- order(o)
        q <- pa <- rep.int( min(n*p/i), n)
        for (j in (n-1):2) {
            ij <- seq_len(n-j+1)
            i2 <- (n-j+2):n
            q1 <- min(j*p[i2]/(2:j))
            q[ij] <- pmin(j*p[ij], q1)
            q[i2] <- q[n-j+1]
            pa <- pmax(pa,q)
        }
        pmax(pa,p)[if(lp < n) ro[1:lp] else ro]
        }
    """
    raise NotImplementedError
#edef

def hochberg(p, n):
    """
    Hochberg correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
    hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin( (n - i + 1L) * p[o] ))[ro]
    }
    """
    lp = len(p)
    i = np.array(range(lp))[::-1] + 1
    o = order(p, decreasing=True)
    ro = order(o)
    return pmin(1, cummin( (n - i + 1) * p[o] ))[ro]
#edef

def BH(p, n):
    """
    Benjamini-Hochberg correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
    BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin( n / i * p[o] ))[ro]
    }
    """
    lp = len(p)
    i = np.array(range(lp))[::-1] + 1
    o = order(p, decreasing=True)
    ro = order(o)
    return pmin(1, cummin( n / i * p[o] ))[ro]
#edef

def BY(p, n):
    """
    Benjamini-Yekutieli correction of p-values
    Inputs:
     - p: An array of p-values
     - n: number of tests to correct for
    Outputs:
     - Corrected q-values
    
    R code:
    BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1L/(1L:n))
    pmin(1, cummin(q * n / i * p[o]))[ro]
    }
    """
    lp = len(p)
    i = np.array(range(lp))[::-1] + 1
    o = order(p, decreasing=True)
    ro = order(o)
    return pmin(1, cummin(q * n / i * p[o]))[ro]
#edef