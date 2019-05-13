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
        for f in self._fields:
            dstr += " %s : %s\n" % (f, self._field_values[f])
        #efor
        dstr += " Values:\n"
        dstr += "  I : %s\n" % ' '.join(self._fields) 
        dstr += '\n'.join(['  %d : %s' % (i, str(r)) for i,r in enumerate(self._values[:10]) ])
        return dstr
    #edef
    
    def __repr__(self):
        return str(self)
    #edef
                
    def __dir__(self):
        """
        To allow tab-completion for the registered objects.
        """
        return object._dir_(self) + list(self._fields) + list(self._all_fields)
    #edef
#eclass

class MappingIndex(object):
    """
    This module allows a mapping between values that have multiple identifiers.
    For example, in BioMart, a gene will have multiple transcripts, and transcripts may have multiple proteins.
    Similarly, multiple proteins may have the same gene identifier
    
    Example usage:
    --------------
    # We define a dataframe which has multiple mappings between different kinds of identifiers, sometimes they are missing.
    # In this case, a mapping between 'small', 'big', and 'und' identifiers:
    df = pd.DataFrame([('a', 'A', None),('a','AA','_a'),('b', 'B', None),('b','BB','_b')], columns=['small', 'big', 'und'])
    
    # Make an index on the 0th column
    mi = MappingIndex(df, key=0)
    
    # You can access the mapping for a specific identifier as if it were a dictionary:
    # It returns a MappingIndexObject
    a_map = mi['a']
    
    You can 
    
    upper = mi['a'].big # -> A
    under = mi['a'].und # -> _a
    other_upper = mi['a'][1].big # -> AA
    
    # You can get a list of all the small, big and und identifiers by prefixing with 'all_':
    mi['a'].all_big
    mi['a'].all_small
    mi['a'].all_und
    
    
    The table used to make the index can be accessed with:
    
    object._table
    """
    

    _slots_ = [ '_idx', '_key', '_file_name', '_tbl', '_names', '_empty_result' ]

    def __init__(self, data, key=0, **kwargs):
        """
        parameters:
        -----------
        data: String|pd.DataFrame. A file to load and index, or a pandas dataframe
        key:  Integer. Which column to use as an index.
        **kwargs: Additional arguments to pd.read_csv
        """
        self._idx = None
        self._file_name = data if isinstance(data, str) else None
        self._tbl = data if isinstance(data, pd.DataFrame) else None
        self._empty_result = None
        self._key = key

        if (self._tbl is None) and (self._file_name is None):
            raise ValueError("You must specify either a filename or a table.")
        elif (self._tbl is None):
            self._tbl = pd.read_csv(self._file_name, **kwargs)
        #fi
        
        self._idx = { k : [] for k in self._tbl[self._tbl.columns[self._key]].unique() }
        for row in self._tbl.values:
            self._idx[row[self._key]].append([ None if pd.isna(v) else v for v in row])
        #efor
        self._idx = { k : MappingIndexObject(self._tbl.columns, v) for (k,v) in self._idx.items()}
        
        self._empty_result = MappingIndexObject(self._tbl.columns, [[ None for c in self._tbl.columns]])
    #edef

    def lookup(self, key):
        """
        lookup: Lookup the value of a key
        Inputs: key :      [str] index value to retrieve
        """
        if key not in self._idx:
            utils.msg.error("Item '%s' not in map." % key)
            return self._empty_result
        #fi

        return self._idx[key]
    #edef
    
    @property
    def _table(self):
        """
        The table used to create the index 
        """
        return self._tbl

    def __contains__(self, key):
        """
        Check if a key is present in the mapping
        """
        return key in self._idx
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
        return self._idx.keys()
    #edef

    def values(self):
        """
        values: Return a list of values of the index
        """
        return self._idx.values()
    #edef

    def __str__(self):
        """
        String representation of the object
        """
        dstr = "Indexed TSV Object\n"
        dstr += " Filename: %s\n" % self._file_name
        dstr += " Indexed on column %s\n" % (str(self._key))
        dstr += " #Rows : %d\n" % self._tbl.shape[0]
        dstr += " #Indexes: %d\n" % len(self._idx)
        return dstr
    #edef
      
    def __repr__(self):
        """
        String representation of the object
        """
        return str(self)
    #edef

#eclass
          
