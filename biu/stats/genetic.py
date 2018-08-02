from .. import utils
from .. import processing

sstats = utils.py.loadExternalModule("scipy.stats")
np     = utils.py.loadExternalModule("numpy")

from collections import namedtuple

def hardyWeinbergEquilibrium(data, alpha=0.001):
    """
      data : Array of 0/1/2 for allele counts per individual
      alpha: alpha cutoff level
      
      output: HardyWeinbergResult namedtuple
    
    Alpha value taken from:
     https://www.nature.com/articles/ng.2756#methods
    TODO: Implement other methods, as defined in
     https://www.jstor.org/stable/pdf/2556115.pdf
    """
    
    data = np.array(data)
    
    resTuple = namedtuple("HardyWeinbergResult", [ 'pvalue', 'method', 'chi2statistic', 'significant' ])

    n   = len(data)
    p   = sum(data) / (sum(data) + sum(data == 0))
    obs = processing.lst.freq(data)
    exp = { 0.0 : np.power(1-p,2) * n,
            1.0 : 2*p*(1-p) * n,
            2.0 : np.power(p,2) * n }

    chi2 = sum([ np.power(obs[c] - exp[c],2)/exp[c] for c in exp ])
    pvalue = 1 - sstats.chi2.cdf(chi2, 1)
    
    print(obs, exp, chi2, pvalue)

    
    return resTuple(pvalue, 'chi2', chi2, pvalue < alpha)
#edef