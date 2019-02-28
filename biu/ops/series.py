from .. import utils

pd = utils.py.loadExternalModule('pandas')
np = utils.py.loadExternalModule('np')

#########################################################################

def cast_str(series):
    """
    Cast a pandas series to numerical type.
    Takes care of nan values.
    
    parameters:
    -----------
    series: A pandas series
    
    Returns:
    --------
    A series with string type
    """
    return series.apply(lambda v: str(v) if str(v).strip() != "" else None)
#edef

#########################################################################

def cast_float(series):
    """
    Cast a pandas series to numerical type.
    Takes care of nan values.
    
    parameters:
    -----------
    series: A pandas series
    
    Returns:
    --------
    A series with numeric type
    """
    def try_cast(value):
        try:
            if not str(value).strip():
                return np.nan
            #fi
            return float(value)
        except Exception as e:
            return np.nan
        #etry
    #edef
    
    return pd.to_numeric(series.apply(try_cast), errors='coerce')
#edef

#########################################################################

def is_categorical(series):
    """
    Detect if a series is a categorical.
    NOT FOOLPROOF. Never is!
    
    Basically, it checks if all non-NA elements in the series are numeric or not.
    If yes, then it is a numerical series.
    If no, then it is a categorical series
    
    parameters:
    -----------
    series: A pandas series
    
    Returns:
    --------
    True if categorical, False otherwise
    """
    if hasattr(series, 'str') or (series.dtype.name == 'object'):
        not_na = series.str.isnumeric()[~pd.isna(series)]
        if all(not_na):
            return False
        else:
            # make it categorical
            return True
        #fi
    else:
        return False
    #fi
#edef
