import collections
import datetime

from .. import utils

rpy2 = utils.py.loadExternalModule('rpy2')
np   = utils.py.loadExternalModule('numpy')

##################################################################################

def dict2ri(D):
    """
    Convert a dictionary to an R ListVector
    """
    if not isinstance(D, dict):
        raise ValueError("Expected dict. Got '%s'." % str(type(D)))
    #fi
    return rpy2.robjects.ListVector(D.items())
#edef

##################################################################################

def ri2dict(D):
    """
    Convert a StrVector to a dictionary
    Note, that this conversion is not the inverse of dict2ri, as R values are always lists...
    Thus, ri2dict(dict2ri({'a': 1})) -> { 'a': [1]}
    """
    
    return dict(zip(D.names, map(list, list(D))))
#edef

##################################################################################

def tuple2ri(T):
    """
    Convert a tuple to an array.
    It is first converted to a numpy array, and then, based on that to an R array
    """
    if not isinstance(T, tuple):
        raise ValueError("Expected tuple. Got '%s'." % str(type(T)))
    #fi
    
    return rpy2.robjects.numpy2ri.numpy2rpy(np.array(T))
#edef

##################################################################################

def none2ri(N):
    """
    Convert a None type to NULL
    """
    return rpy2.robjects.NULL
#edef

##################################################################################

def datetime2ri(D):
    """
    Convert a datetime object to an R object
    """
    import rpy2.robjects.packages as rpackages
    base = rpackages.importr("base") 
    rdate = base.as_POSIXlt(D.strftime("%Y-%m-%d %H:%M:%S"), format="%Y-%m-%d %H:%M:%S")
    return rdate
#edef

##################################################################################

def converter():
    """
    Return an rpy2 converter that automatically converts several formats
    
    Automatically added:
    * pandas converted
    * numpy objects
    * dict
    * tuple
    """
    from rpy2.robjects import numpy2ri
    from rpy2.robjects import pandas2ri
    
    my_converter  = rpy2.robjects.conversion.Converter('BIU converter')
    
    my_converter.rpy2py.register(rpy2.rinterface.ListSexpVector, lambda x: x)
    
    my_converter.py2rpy.register(dict, dict2ri)
    my_converter.rpy2py.register(rpy2.robjects.ListVector, ri2dict)
    
    my_converter.py2rpy.register(tuple, tuple2ri)
    my_converter.py2rpy.register(type(None), none2ri)
    my_converter.py2rpy.register(datetime.datetime, datetime2ri)
    
    my_converter += rpy2.robjects.default_converter
    my_converter += numpy2ri.converter
    my_converter += pandas2ri.converter

    return my_converter
#edef

##################################################################################