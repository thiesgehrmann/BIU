from .. import utils

np = utils.py.loadExternalModule('numpy')
sstats = utils.py.loadExternalModule('scipy', 'stats')

#################################################################################

def q_(x, dist=sstats.norm, *pargs, **kwargs):
    """
    A wrapper to make the qnorm, qchi qbinom etc in R
    """
    return dist.ppf(x, *pargs, **kwargs)
#edef

#################################################################################

def p_(x, dist=sstats.norm, *pargs, **kwargs):
    """
    A wrapper to make the pnorm, pchi pbinom etc in R
    """
    return dist.cdf(x, *pargs, **kwargs)
#edef

#################################################################################

def r_(*pargs, dist='normal', **kwargs):
    """
    A wrapper to make the rnorm, rchi rbinom etc in R
    """
    return getattr(np.random, dist)(*pargs, **kwargs)
#edef

#################################################################################