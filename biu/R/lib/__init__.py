from .. import utils
from .. import settings


def modules():
    """
    Provides an index of the R code library in the biu package
    """
    inst_loc = settings.biuLocation() + '/biu/R/src'
    
    D = { "analysis.rnaseq.diffex" : '%s/analysis.rnaseq.diffex.R' }
    
    return D
#edef