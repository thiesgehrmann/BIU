from .. import settings
from .. import utils

import os
import pathlib

class DataObjects(object):
    """
    An object that manages the lazy loading of objects and their file dependencies
    
    Example usage:
    --------------
    
    import biu
    import pandas as pd

    file = biu.utils.Acquire().curl('https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data')
    data_objs = DataObjects(where='./', redo=False)
    data_objs.add_file("iris.tsv", file)
    data_objs.register("iris", ["iris.tsv"],
                       lambda x: pd.read_csv(x["iris.tsv"], index_col=False, names=['a','b','c','d','class']),
                       docstring="An Pandas DataFrame of the IRIS data")
                       
    data_objs.iris
    data_objs.iris['e'] = data_objs.iris.a
    """
    
    __slots__ = [ 'registered', 'files', 'local_files', 'loaded', 'where',
                  'download_where', 'redo', 'acquired_files' ]

    
    def __init__(self, where, download_where, redo=False, local_files=None):
        """
        Initialize a DataObject object
        Parameters:
        -----------
        where:          String. Where the final data should be stored
        download_where: String. Where the data should be downloaded
        redo:           Boolean. Re-download the data, or not
        local_files:    Dict: A dictionary of (name:path) for files
        add_property:   Python class. If specified, objects will be added as properties to this class.
        """
        
        self.where          = where
        self.download_where = download_where
        self.redo           = redo
        self.files          = {}
        self.local_files    = {}
        self.acquired_files = []
        self.registered     = {}
        self.loaded         = {}
        
        if local_files is not None:
            for lf_name, lf_path in local_files.items():
                self.local_files[lf_name] = utils.Acquire2(lf_path)
                self.files[lf_name]       = self.local_files[lf_name]
            #efor
        #fi
        
    #edef
    
    def add_file(self, name, acquire_object, redo=None, finalize=True):
        """
        Add a file acquire object.
        Parameters:
        -----------
        
        name:           String. The name of the file
        acquire_object: biu.utils.Acquire2. The file acquisition pipeline
        redo:           boolean. Set the redo status of this file, regardless of the total data object (Default is None)
        finalize:       boolean|String
                        if bool and true: at the end of the file acquisition, copy the file to "self.where/name".
                        if string: at the end of the file acquisition, copy the file to "self.where/finalize".
                        
        """
        if name in self.local_files:
            utils.msg.dbm("The file '%s' has been specified by a localCopy, I will not add the new file specification." % name)
        else:
            redo = self.redo if redo is None else redo
            
            if not isinstance(acquire_object, utils.Acquire2):
                raise ValueError("Expected a biu.Acquire2 object.")
            #fi
            
            acquire_object = acquire_object.set_redo(redo).set_where(self.download_where)
            
            if isinstance(finalize, str):
                acquire_object = acquire_object.finalize("%s/%s" % (self.where, finalize))
            elif isinstance(finalize, bool) and finalize:
                acquire_object = acquire_object.finalize("%s/%s" % (self.where, name))
            else:
                acquire_object = acquire_object.copy()
            #fi
            
            self.files[name] = acquire_object
        #fi
    #edef
    
    def get_file(self, name):
        """
        Return the acquire file object.
        """
        if name not in self.files:
            raise ValueError("No such file is known '%s'." % name)
        #fi
        return self.files[name]
    #edef
    
    def _acquire_files(self, files):
        """
        Perform the acquisition pipeline for the files specified.
        parameters:
        -----------
        files: List of strings
        """
        for f in files:
            if f in self.acquired_files:
                continue
            #fi

            self.files[f].acquire()
        #efor
    #edef
    
    def register(self, name, required_files, load_func):
        """
        Register an object
        
        parameters:
        -----------
        name:           String. python-identifier valid name of the object
        required_files: list. List of files (added by add_files) that need to exist before loading data
        load_func:      Python function. Takes as input dictionary of paths. Must output some object
        docstring:      String.
        
        NOTE: One or more of your arguments will most likely be the file of one of your files.
        If you have not specified a finalized location of these files, you must pass a placeholder parameter using file_path(file_name) such as this:
        
        Example:
        --------
        
        file = biu.utils.Acquire2().curl('https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data')
        data_objs = DataObjects(where='./', redo=False)
        data_objs.add_file("iris.tsv", file)
        data_objs.register("iris", ["iris.tsv"],
                           lambda x: pd.read_csv(x["iris.tsv"], index_col=False, names=['a','b','c','d','class']),
                           docstring="An Pandas DataFrame of the IRIS data")
        """
        if not name.isidentifier():
            raise ValueError("The name '%s' is not a valid object name. Must be a valid python identifier." % name)
        #fi
        
        if name in self.__slots__:
            raise ValueError("Sorry, '%s' is not a valid object name. It is already used by this class." % name)
        #fi
        
        for file in required_files:
            if file not in self.files:
                raise ValueError("No such file is known: '%s'." % file)
            #fi
        #efor
        
        if not hasattr(load_func, '__call__'):
            raise ValueError('The provided load_func is not callable.')
        #fi
        
        self.registered[name] = (load_func, required_files)
        
        getter = lambda x: self.get(name)
        setter = lambda x, value: self.set(name, value)
        prop   = property(fget=getter, fset=setter)
        setattr(self.__class__, name, prop)
    #edef
    
    def load(self, name):
        """
        Loads or re-loads an object
        """
        if name not in self.registered:
            raise AttributeError("No such object is registered: '%s'." % name)
        #fi
        
        load_func, required_files = self.registered[name]
        self._acquire_files(required_files)
        
        self.loaded[name] = load_func({ file_name : self.files[file_name].path for file_name in required_files })
    #edef
    
    def isloaded(self, name):
        if name not in self.registered:
            raise NameError("No such object is registered: '%s'." % name)
        #fi
        return name in self.loaded
    #edef
    
    def get(self, name):
        """
        Return an object. Load if not loaded. Raises ValueError if object is not known
        """
        if name not in self.loaded:
            self.load(name)
        #fi
        
        return self.loaded[name]
    #edef
    
    def __getitem__(self, name):
        """
        Return an object. Load if not loaded. Raises ValueError if object is not known
        """
        return self.get(name)
    #edef
    
    def set(self, name, value):
        """
        Change the value of a loaded object (Useful if you are doing data management within the dataset)
        """
        if (name not in self.loaded) and (name not in self.registered):
            raise AttributeError("No such object is known: '%s'." % name)
        #fi
        self.loaded[name] = value
    #edef
    
    def __setitem__(self, name, value):
        """
        Change the value of a loaded object (Useful if you are doing data management within the dataset)
        """
        self.set(name, value)
    #edef

    def __contains__(self, name):
        """
        Check if an object is registered.
        """
        return name in self.registered
    #edef
    
#eclass

##############################################################

def Dataset2_test():
    """
    A wrapper function that makes a new version of the Dataset2 class when it is run.
    """
    
    class Dataset2_unique_class(Dataset2):
        def __init__(self, *pargs, **kwargs):
            super(Dataset2_unique_class, self).__init__(*pargs, **kwargs)
    #eclass
    
    return Dataset2_unique_class
#edef

class Dataset2(object):
    """
    A BIU Dataset Object (Version 2).
    Datasets consist of files and objects.
    The Dataset object manages the dynamic acquisition of files and the dynamic loading of files into objects.
    For example, the IRIS dataset consists of a csv file, which is loaded into a dataframe object.
    
    Datasets will also provide additional methods to specifically access the relevant dataset
    """
    
    __slots__ = [ '__where', '__download_where', '__redo', '_obj', '_str_funcs' ]
    
    def __init__(self, dataset_identifier,
                 where=None,
                 download_where=None,
                 redo=False, local_files=None):
        """
        Initialize a Dataset object:
        Parameters:
        -----------
        dataset_identifier: String. A name for the dataset. Typically 'dataset_name/version'
        where:              String. Where the final data should be stored
        download_where:     String. Where the data should be downloaded
        redo:               Boolean. Re-download the data, or not
        local_files:        Dict: A dictionary of (name:path) for files
        """
        
        where = where if where is not None else settings.getDataDir()
        download_where = download_where if download_where is not None else settings.getDownloadDir()
        
        class _data_class(DataObjects):
            def __init__(self, *pargs, **kwargs):
                super(_data_class, self).__init__(*pargs, **kwargs)
        #eclass
        
        for attr in object.__getattribute__(self, '__slots__'):
            object.__setattr__(self, attr, None)
        #efor
        
        where          = os.path.abspath(os.path.expanduser(where) + '/' + dataset_identifier)
        download_where = os.path.abspath(os.path.expanduser(download_where))
        
        object.__setattr__(self, '__where', where)
        object.__setattr__(self, '__download_where', download_where)
        object.__setattr__(self, '__redo', redo)
        object.__setattr__(self, '_obj', _data_class(where, download_where, redo, local_files))
        object.__setattr__(self, '_str_funcs', [])
    #edef
    
    def _add_str_func(self, fun):
        """
        Add a function that is evaluated and printed when a string representation is made.

        Parameters:
        -----------
        fun: Function. Dataset2 object -> string.
             e.g. lambda x: '\n'.join(list(x._obj.files.keys()))
        """
        
        if not hasattr(fun, '__call__'):
            raise ValueError("The item provided must be callable.")
        #fi
        
        self._str_funcs.append(fun)
    #edef
        
    
    def __str__(self):
        """
        Prepare a string representation of the class
        """
        dstr  = "%s object\n" % self.__class__.__name__
        
        dstr += "Where: %s\n" % object.__getattribute__(self, '__where')
        dstr += "Downloaded data in: %s\n" % object.__getattribute__(self, '__download_where')
        dstr += "Redo: %s\n" % str(object.__getattribute__(self, '__redo'))

        for f in self._str_funcs:
            fstr = f(self)
            if fstr is None:
                continue
            #fi
            for line in fstr.split('\n'):
                dstr += ' ' + line + '\n'
            #efor
        #efor

        dstr += ' Objects:\n'
        for oname in self._obj.registered:
            if oname[0] == '_':
                continue
            #fi
            loaded = oname in self._obj.loaded
            dstr += '  * [%s] %s\n' % (('X' if loaded else ' '), oname)
        #efor

        dstr += " Files:\n"
        for what in self._obj.files:
            loc = self._obj.get_file(what).path
            if os.path.islink(loc):
                dstr += "  * [%s] %s : %s -> %s\n" % ('S' if self._obj.get_file(what).exists else ' ', what, loc, pathlib.Path(loc).resolve())
            else:
                dstr += "  * [%s] %s : %s\n" % ('X' if self._obj.get_file(what).exists else ' ', what, loc)
            #fi
        #efor
        return dstr
    #edef
    
    def __repr__(self):
        """
        Prepare a string representation of the class
        """
        return str(self)
    #edef
    
    def __getattr__(self, name):
        """After regular attribute access, try looking up the name within the registered objects
        This allows simpler access to objects for interactive use.
        """

        # Note: obj.x will always call obj.__getattribute__('x') prior to
        # calling obj.__getattr__('x').
        
        try:
            object.__getattribute__(self, name)
        except AttributeError:
            pass
        #etry
        
        return object.__getattribute__(self, '_obj')[name]

    #edef

    def __setattr__(self, name, value):
        """After regular attribute access, try setting the name
        This allows simpler access to objects for interactive use.
        """

        # first try regular attribute access via __getattribute__, so that
        # e.g. ``obj.x`` and ``obj.x = 4`` will always reference/modify
        # the same attribute.

        try:
            object.__getattribute__(self, name)
            return object.__setattr__(self, name, value)
        except AttributeError:
            pass
        #etry
        
        object.__getattribute__(self, '_obj')[name] = value
    #edef
    
    def __dir__(self):
        """
        To allow tab-completion for the registered objects.
        """
        return object.__dir__(self) + list(self._obj.registered.keys())
    #edef

#eclass