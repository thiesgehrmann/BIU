from .. import utils
from .. import stats
from .. import ops

pd = utils.py.loadExternalModule("pandas")

class GAF(object):
    """
    A GO Annotation Format handler.
    """

    # GO annotation File: http://geneontology.org/page/go-annotation-file-gaf-format-21
  
    _fieldNames = [ "db", "db_o_id", "db_o_symbol", "qualifier", "go_id", "db_ref", "evidence", "w_o_f", "aspect", "db_o_name", "db_o_synonym", "db_o_type", "date", "assigned_by", "annot_ext", "gene_product_form_id" ]
  
    _entries = None
    _fileName = None
    _annotIndex = None
    _objectIndex = None
  
    def __init__(self, file_name, **kwargs):
        """
        Initialize the GAF object.
        
        parameters:
        -----------
        file_name: String. A path to the file
        **kwargs: Additional arguments to pd.read_csv
        """
        self._entries = pd.read_csv(file_name, index_col=False, names=self._fieldNames, **kwargs).drop_duplicates(('go_id', 'db_o_id'))
        self._fileName = file_name
        self._annotIndex = {}
        self._objectIndex = {}
  
        for i, row in self._entries.iterrows():
            goID = row.go_id
            objectID = row.db_o_id
    
            if goID not in self._annotIndex:
                self._annotIndex[goID] = []
            #fi
            if objectID not in self._objectIndex:
                self._objectIndex[objectID] = []
            #fi
    
            self._annotIndex[goID].append( (objectID, i) )
            self._objectIndex[objectID].append( (goID, i) )
        #efor
    #edef
  
    @property
    def annotations(self):
      """
      Return the list of all annotations.
      """
      return list(self._annotIndex.keys())
    #edef
  
    @property
    def objects(self):
      """
      Return a list of all annotated objects.
      """
      return list(self._objectIndex.keys())
    #edefd
  
    @utils.decorators.deprecated("getAnnots is deprecated. Use get_annots instead.")
    def getAnnots(self, objectID):
        """
        Return the annotations for a given object.
        Returns empty list if object is unknown.
        """
        return self.get_annots(objectID)
    #edef
    
    def get_annots(self, object_id):
        """
        Return the annotations for a given object.
        Returns empty list if object is unknown.
        """
        if object_id not in self._objectIndex:
            return []
        else:
            return [ a[0] for a in self._objectIndex[object_id] ]
        #fi
    #edef
    
    @utils.decorators.deprecated("getAnnoted is deprecated. Use get_annots instead.")
    def getAnnotated(self, annotID):
        """
        Return the annotations for a given object.
        Returns empty list if annotation is unknown
        """
        return self.get_annotated(annotID)
    #edef
  
    def get_annotated(self, annot_id):
        """
        Return the annotations for a given object.
        Returns empty list if annotation is unknown
        """
        if annot_id not in self._annotIndex:
            return []
        else:
            return [ a[0] for a in self._annotIndex[annot_id] ]
        #fi
    #edef
  
    def enrich(self, your_set, pathway=None, method=None, **kwargs):
        """
        Enrich: Check enrichment of GO pathways in a given set
  
        Inputs:
          - yourSet: List of Uniprot protein IDs to test
          - pathway: List of GO terms (or single term) to test (Defaults to all terms that your objectIDs are present in)
          - method: Type of multiple testing correction procedure to use
          - **kwargs: Additional arguments for multple testing procedure
  
        Outputs:
         - df : Pandas Data Frame of test results
        """
        if pathway is None:
            pathway = list(set([ p for prot in your_set for p in self.get_annots(prot) ]))
        elif isinstance(pathway, str):
            pathway = [ pathway ]
        #fi
  
        R = []
        B = self.objects
        for p in pathway:
            pathway_genes = self.get_annotated(p)
            res = stats.enrichment.set_enrichment(your_set, pathway_genes, B, abcd_values=abcd_values)
            R.append((p, res.method, res.c2statistic, res.oddsratio, res.pvalue, res.table, res.table_values))
        #efor
  
        df = pd.DataFrame(R, columns=['pathway', 'method', 'c2statistic', 'oddsratio', 'p', 'table', 'table_values'])
        if method is not None:
            df['q'] = stats.p_adjust(df.p.values, method, **kwargs)
        #fi
  
        return df
    #edef
  
    def summary(self, your_sets):
        """
        For a set of proteins, count the occurrence of each pathway.
  
        parameters:
        -----------
        your_sets: List[String] | List[List[String]] | dict[string->List[String]]
            A list of uniprot protein identifiers, or a list/dict of lists of uniprot identifiers
  
        Returns:
        Pandas Dataframe of pathway counts
          
        example usage:
          
        """
  
        if not isinstance(your_sets, dict):
            if isinstance(your_sets[0], str):
                your_sets = { "count" : your_sets }
            else:
                your_sets = { "count_%d" % (i+1) : s for (i,s) in enumerate(your_sets) }
            #fi
        #fi
          
        tot_set = set(self.objects)
        your_sets = { k : list(set(your_sets[k]) & tot_set) for k in your_sets }
          
        p_counts = []
          
        for this_k in your_sets:
            this_set = your_sets[this_k]
            pathways = [ pway for prot in this_set for pway in self.get_annots(prot) ]
            pathways = [ (this_k, k, v) for (k,v) in  ops.lst.freq(pathways).items() ]
            p_counts.extend(pathways)
        #efor
        p_counts = pd.DataFrame(p_counts, columns=[ 'set', 'pathway', 'count'])
        p_counts = pd.pivot_table(index='pathway', columns='set', values='count', data=p_counts)
        p_counts = p_counts.fillna(0)
        for this_k in your_sets:
            if this_k not in p_counts.columns:
                p_counts[this_k] = [0] * p_counts.shape[0]
            #fi
            p_counts[this_k] = p_counts[this_k].astype(int)
        #efor
        p_counts['total_count'] = [ len(self.get_annotated(p)) for p in p_counts.index ]
        #p_counts['description'] = [ self.pathway[p].description for p in p_counts.index ]
  
  
        return p_counts
    #edef
      
    def annotate(self, your_sets):
        """
        For a set of proteins, Annotate them
  
        parameters:
        -----------
        your_sets: List[String] | List[List[String]] | dict[string->List[String]]
            A list of uniprot protein identifiers, or a list/dict of lists of uniprot identifiers
  
        Returns:
        Pandas Dataframe of pathway counts
          
        example usage:
          
        """
  
        if not isinstance(your_sets, dict):
            if isinstance(your_sets[0], str):
                your_sets = { "annot" : your_sets }
            else:
                your_sets = { "annot_%d" % (i+1) : s for (i,s) in enumerate(your_sets) }
            #fi
        #fi
          
        tot_set = set(self.objects)
        your_sets = { k : set(your_sets[k]) & tot_set for k in your_sets }
          
        p_annot = []
          
        pathways = list(set([ pway for p in ops.lst.flatten(your_sets.values()) for pway in self.get_annots(p) ]))
          
        for pway in pathways:
            pway_set = set(self.get_annotated(pway))
            counts = [ your_sets[k] & pway_set for k in your_sets ]
            p_annot.append([pway] + counts)
        #efor
          
          
        p_annot = pd.DataFrame(p_annot, columns=[ 'pathway' ] + list(your_sets.keys()))
  
        return p_annot
    #edef
  
    def __str__(self):
        """
        String representation of object.
        """
        dstr  = "GAF (GO Annotation File) Object\n"
        dstr += " Where: %s\n" % self._fileName
        dstr += " # Annotations : %d\n" % self._entries.shape[0]
        dstr += " # Objects     : %d\n" % len(self._objectIndex.keys())
        dstr += " # GO terms    : %d\n" % len(self._annotIndex.keys())
        return dstr
    #edef
    
    def __repr__(self):
        """
        String representation of object.
        """
        return str(self)
    #edef

#eclass
