from ..structures import Dataset2
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule('pandas')

class MultiMappingIndex(Dataset2):
    """
    For a table with multiple indexes, construct a multimapping index.
    Same interface as the BioMart interface, but for any text file.
    
    The internal construction is a Dataset2, but this is hidden.
    
    Usage:
    
    MultiMappingIndex(object, *pargs, **kwargs)
    
    parameters:
    -----------
    object : String|pandas.DataFrame
        The location of a file, OR a pre-specified pandas dataframe
        *pargs, **kwargs: arguments to pd.read_csv

    """
    
    def __init__(self, obj, keys=None, *pargs, **kwargs):
        super(MultiMappingIndex, self).__init__("MultiMappingIndex")
        
        if isinstance(obj, str):
            obj = pd.read_csv(obj, *pargs, **kwargs)
        #fi
        
        if not isinstance(obj, pd.DataFrame):
            raise Exception('Object is not a file or a pandas dataframe')
        #fi
        
        import re
        def clean(s):
            s = str(s)
            # Remove invalid characters
            s = re.sub('[^0-9a-zA-Z_]', '', s)
            
            if s[0] in '1234567890':
                s = 'c' + s
            #fi
            
            return s
        #edef
        
        obj = obj.copy()
        obj = obj.rename(columns={c: clean(c) for c in obj.columns})
        
        object.__setattr__(self, '_table', obj)
        
        if keys is None:
            keys = obj.columns
        #fi
        
        for idx, colname in enumerate(obj.columns):
            if colname in keys:
                self._obj.register(colname, [], lambda f, idx=idx: formats.MappingIndex(self._table, idx))
            #fi
        #efor
    #edef
    
    
    
    def __str__(self):
        """
        Prepare a string representation of the class
        """
        dstr  = "%s object\n" % self.__class__.__name__

        dstr += ' Indexes:\n'
        for oname in self._obj.registered:
            if oname[0] == '_':
                continue
            #fi
            loaded = oname in self._obj.loaded
            dstr += '  * [%s] %s\n' % (('X' if loaded else ' '), oname)
        #efor

        return dstr
    #edef
#eclass