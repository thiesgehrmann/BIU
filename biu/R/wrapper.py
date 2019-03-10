from .. import utils
rpy2 = utils.py.loadExternalModule('rpy2')
np   = utils.py.loadExternalModule('numpy')
pd   = utils.py.loadExternalModule('pandas')



class R(object):
    """
    A wrapper for rpy2, which somewhat mimics the ipython magic functions.
    Basically, it handles the automatic conversion of some python objects to R objects.
    Further, it allows you to automatically push python objects, call code and get converted objects back to python.
    
    Example usage:
    --------------
    
    r = biu.R()
    x = pd.DataFrame([[1,2,3],[4,5,6]])
    r.push(x=x)
    r('y = x * 2')
    y = r.get('y')
    
    Or, altogether:
    ---------------
    y = r('y=x*2', push=dict(x=x), get='y')
    
    Doing a lot at the same time:
    -----------------------------
    
    y, z = r('''
        y = x * 2
        z = x + 2
        ''', push=dict(x=x), get=['y', 'z'])
    
    
    """
    _tconverter          = None
    _converter           = None
    
    def __init__(self):
        """
        Initialize the rpy2 wrapper
        """

        from rpy2.robjects.conversion import converter as template_converter

        from rpy2.robjects import numpy2ri
        from rpy2.robjects import pandas2ri
        template_converter += numpy2ri.converter
        template_converter += pandas2ri.converter

        self._tconverter = template_converter
        self._converter  = rpy2.robjects.conversion.Converter('BIU converter', template=template_converter)
    #edef
    
    def add_converter(self, converter):
        """
        Add a converter to the object, if there is one missing.
        """
        self._tconverter += converter
        self._converter  = rpy2.robjects.conversion.Converter('BIU converter', template=template_converter)
    #edef
        
    def push(self, **kwargs):
        """
        Push values to R, based on the current converter
        
        parameters:
        -----------
        kwargs: Dictionary of values
        
        Example usage:
        --------------
        
        r.push(x=10, y='pool', ages=[10, 50, 100])
        """
        
        if kwargs is None:
            return None
        #fi
        
        for (k,v) in kwargs.items():
            with rpy2.robjects.conversion.localconverter(self._converter) as cv:
                rpy2.robjects.r.assign(k, v)
            #ewith
        #efor
    #edef
        
    def get(self, name, *pargs):
        """
        Get a value from R, based on the current converter
        
        parameters:
        -----------
        name: return this variable from the R instance
        *pargs, if specified, return a tuple of name + those in pargs
        
        returns:
        --------
        Either a converted R object, or
        if pargs is specified, then a tuple of values
        """
        with rpy2.robjects.conversion.localconverter(self._converter):
            if len(pargs)  == 0:
                return rpy2.robjects.globalenv.find(name)
            else:
                return [ rpy2.robjects.globalenv.find(n) for n in ([name] + list(pargs)) ]
            #fi
        #ewith
            
    #edef
    
    def exec(self, cmd, push=None, get=None):
        """
        Call R code, pushing values, and returning values if necessary
        
        parameters:
        -----------
        
        cmd: The R code you want to execute
        push: Dictionary of name:value pairs that you want to introduce to R session
        get: List of R object values that you want to get back
        
        returns:
        A value, or tuple of the get parameters you specified.
        If you specified a single get
        
        """
        if push is None:
            push = {}
        #fi
        
        self.push(**push)
        
        rpy2.robjects.r(cmd)
        
        if get is None:
            return None
        elif isinstance(get, str):
            return self.get(get)
        else:
            return self.get(*get)
        #fi
    #edef
    
    def __call__(self, cmd, push=None, get=None):
        """
        Call R code, pushing values, and returning values if necessary
        
        parameters:
        -----------
        
        cmd: The R code you want to execute
        push: Dictionary of name:value pairs that you want to introduce to R session
        get: List of R object values that you want to get back
        
        returns:
        A value, or tuple of the get parameters you specified.
        If you specified a single get
        
        """
        return self.exec(cmd, push, get)
    #edef
#eclass