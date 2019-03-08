"""
This module provides many functions available in R. This is to make it easier to convert code between R and python. These functions 

At some point, I hope it will contain versions of all the following methods:
http://www.sr.bham.ac.uk/~ajrs/R/r-function_list.html
"""

from ..ops.array import order, pmin, pmax, cummin, cummax

from ..stats.p_adjust import p_adjust

from ..stats.regression import lowess

from .. import settings
from .. import utils

rpy2 = utils.py.loadExternalModule('rpy2')

@property
def modules():
    """
    Provides an index of all R code files in the biu package.
    
    Provided as an rpy2 ListVector
    """
    inst_loc = settings.biuLocation() + '/biu/R/src'
    
    D = { "analysis.rnaseq.diffex" : '%s/analysis.rnaseq.diffex.R' }
    
    return rpy2.robjects.ListVector({ m : D[m] % inst_loc for m in D })
#edef

from .wrapper import R