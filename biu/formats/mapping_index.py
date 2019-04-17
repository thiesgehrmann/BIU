from .. import utils

pd = utils.py.loadExternalModule("pandas")

class MappingIndexObject(object):
    def __init__(self, fields, values):
        """
        Make a MappingIndexObject
        
        parameters:
        -----------
        fields: List[String]. the named fields to use
        values: 
        """
        
        nfields = len(fields)
        
        if not all([len(x) == nfields for x in values]):
            raise ValueError("Incorrect number of values for specified fields.")
        #fi
        
        self._fields = fields
        self._all_fields = [ 'all_%s' % f for f in fields ]
        self._values = values
        
        fv = [ [ v for v in set(f) if v is not None ] for f in zip(*values) ]
        self._all_field_values = { f : v for (f,v) in zip(self._all_fields, fv)}
        
        self._field_values = { f : v[0] if len(v) > 0 else None for (f,v) in zip(fields, fv) }
    #edef
    
    def __getitem__(self, index):
        return MappingIndexObject(self._fields, [self._values[index]])
    #edef
    
    def lookup(self, attr):
        """
        Lookup the value of a field.
        parameters:
        -----------
        attr: String. The name of the field
        """
        if attr in self._fields:
            return self._field_values[attr]
        elif attr in self._all_fields:
            return self._all_field_values[attr]
        else:
            raise AttributeError("Unknown field '%s'." % attr)
        #fi
    #edef
    
    def __getattr__(self, attr):
        if attr in self._fields:
            return self._field_values[attr]
        elif attr in self._all_fields:
            return self._all_field_values[attr]
        else:
            raise AttributeError("Unknown attribute '%s'." % attr)
        #fi
    #edef
    
    def __str__(self):
        dstr  = "MappingIndexObject\n"
        for f in self._fields[:10]:
            dstr += " %s : %s\n" % (f, self._field_values[f])
        #efor
        dstr += " Values:\n"
        dstr += "  I : %s\n" % ' '.join(self._fields) 
        dstr += '\n'.join(['  %d : %s' % (i, str(self._values[i])) for i in range(len(self._values)) ])
        return dstr
    #edef
    
    def __repr__(self):
        return str(self)
    #edef
                
    def __dir__(self):
        """
        To allow tab-completion for the registered objects.
        """
        return object.__dir__(self) + list(self._fields) + list(self._all_fields)
    #edef
#eclass

class MappingIndex(object):
    """
    This module allows a mapping between values that have multiple identifiers.
    For example, in BioMart, a gene will have multiple transcripts, and transcripts may have multiple proteins.
    Similarly, multiple proteins may have the same gene identifier
    """
    

    __slots__ = [ '__idx', '__key', '__file_name', '__tbl', '__names', '__empty_result' ]

    def __init__(self, data, key=0, **kwargs):
        """
        parameters:
        -----------
        data: String|pd.DataFrame. A file to load and index, or a pandas dataframe
        key:  Integer. Which column to use as an index.
        **kwargs: Additional arguments to pd.read_csv
        """
        self.__idx = None
        self.__file_name = data if isinstance(data, str) else None
        self.__tbl = data if isinstance(data, pd.DataFrame) else None
        self.__empty_result = None
        self.__key = key

        if (self.__tbl is None) and (self.__file_name is None):
            raise ValueError("You must specify either a filename or a table.")
        elif (self.__tbl is None):
            self.__tbl = pd.read_csv(self.__file_name, **kwargs)
        #fi

        
        
        self.__idx = { k : [] for k in self.__tbl[self.__tbl.columns[self.__key]].unique() }
        for row in self.__tbl.values:
            self.__idx[row[self.__key]].append([ None if pd.isna(v) else v for v in row])
        #efor
        self.__idx = { k : MappingIndexObject(self.__tbl.columns, v) for (k,v) in self.__idx.items()}
        
        self.__empty_result = MappingIndexObject(self.__tbl.columns, [[ None for c in self.__tbl.columns]])
    #edef

    def lookup(self, key):
        """
        lookup: Lookup the value of a key
        Inputs: key :      [str] index value to retrieve
        """
        if key not in self.__idx:
            utils.msg.error("Item '%s' not in map." % key)
            return self.__empty_result
        #fi

        return self.__idx[key]
    #edef

    def __contains__(self, key):
        """
        Check if a key is present in the mapping
        """
        return key in self.__idx
    #edef

    def __getitem__(self, key):
        """
        
        Get a specific mapping"""
        return self.lookup(key)
    #edef


    def keys(self):
        """
        keys: Return a list of keys of the index
        """
        return self.__idx.keys()
    #edef

    def values(self):
        """
        values: Return a list of values of the index
        """
        return self.__idx.values()
    #edef

    def __str__(self):
        """
        String representation of the object
        """
        dstr = "Indexed TSV Object\n"
        dstr += " Filename: %s\n" % self.__file_name
        dstr += " Indexed on column %s\n" % (str(self.__key))
        dstr += " #Rows : %d\n" % self.__tbl.shape[0]
        dstr += " #Indexes: %d\n" % len(self.__idx)
        return dstr
    #edef
      
    def __repr__(self):
        """
        String representation of the object
        """
        return str(self)
    #edef

#eclass
          
