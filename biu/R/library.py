import os

from .. import settings

###############################################################################################

def library():
    """
    Provides an index of all R code files in the biu package.
    
    These can be used in the rpy2 interface.
    
    example usage:
    --------------
    
    r = biu.R()
    r('names(biu)')
    r('source(biu$limma.utils)')
    """
    lib_loc = settings.biuLocation + '/biu/R/lib'
    
    D = { '.'.join(f.split('.')[:-1]) : '%s/%s' % (lib_loc, f) for f in os.listdir(lib_loc)
            if (f.split('.')[-1] == 'R') and (f.split('.')[0] != '') }
    
    return D
#edef

###############################################################################################

def data():
    pass
#edef