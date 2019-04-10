from ..structures import Dataset2
from .. import utils
from .. import stats
from .. import ops

pd = utils.py.loadExternalModule('pandas')

class Reactome(Dataset2):
    """
    A class to manage the Reactome Dataset.
    
    example usage:
    
    r = biu.db.Reactome('Homo sapiens')
    proteins = [ prot.id for prot in r.pathway['R-HSA-8956319'].proteins ]
    r.enrich(proteins)
    
    """
    
    _organisms = None
    
    def __init__(self, organisms=None, *pargs, **kwargs):
        """
        Initialize a Reactome instance.
        
        parameters:
        -----------
        organisms: String|list[Strings]. Specify which organisms should be indexed in this instance
        *pargs, **kwargs. See documentation of biu.structures.Dataset2 for details
        """
        super(Reactome, self).__init__("Reactome", *pargs, **kwargs)
        
        if isinstance(organisms, str):
            organisms = [ organisms ]
        #fi
        
        self._organisms = organisms
        
        self._obj.add_file('pathway_uniprot.tsv', utils.Acquire2().curl('https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt'))
        self._obj.add_file('pathway_chemebi.tsv', utils.Acquire2().curl('https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt'))
        self._obj.add_file('pathway_names.tsv', utils.Acquire2().curl('https://reactome.org/download/current/ReactomePathways.txt'))
        self._obj.add_file('pathway_hierarchy.tsv', utils.Acquire2().curl('https://reactome.org/download/current/ReactomePathwaysRelation.txt'))
        
        self._obj.register("pathway_protein_map", ['pathway_uniprot.tsv'],
                           lambda f: pd.read_csv(f['pathway_uniprot.tsv'], sep='\t',
                                                 names=['uniprot', 'pathway', 'url', 'name', 'evidence', 'organism']))
        self._obj.register("pathway_metabolite_map", ['pathway_chemebi.tsv'],
                           lambda f: pd.read_csv(f['pathway_chemebi.tsv'], sep='\t',
                                                 names=['chemebi', 'pathway', 'url', 'name', 'evidence', 'organism']))
        self._obj.register("pathway_names", ['pathway_names.tsv'],
                           lambda f: pd.read_csv(f['pathway_names.tsv'], sep='\t',
                                                 names=['pathway', 'description', 'organism']))
        self._obj.register("pathway_hierarchy", ['pathway_hierarchy.tsv'],
                           lambda f: pd.read_csv(f['pathway_hierarchy.tsv'], sep='\t',
                                                 names=['parent', 'child']))
        
        def make_index():
            """
            Load the reactome into a useful structure
            """

            from collections import namedtuple

            pathways = {}
            proteins = None
            metabolites = None
            
            pathway_names           = self.pathway_names if self._organisms is None else self.pathway_names[self.pathway_names.organism.isin(self._organisms)]
            pathway_protein_map     = self.pathway_protein_map if self._organisms is None else self.pathway_protein_map[self.pathway_protein_map.organism.isin(self._organisms)]
            pathway_metabolite_map = self.pathway_metabolite_map if self._organisms is None else self.pathway_metabolite_map[self.pathway_metabolite_map.organism.isin(self._organisms)]

            print('Adding pathways')
            for i, p in pathway_names.iterrows():
                pathways[p.pathway] = self.ReactomePathway(p.pathway, p.description, p.organism)
            #efor

            print('Adding proteins')
            proteins = { p : self.ReactomeNode(p, 'protein') for p in pathway_protein_map.uniprot }
            for i, a in pathway_protein_map.iterrows():
                if a.pathway not in pathways:
                    pathways[a.pathway] = self.ReactomePathway(a.pathway, 'Unknown pathway', a.organism)
                #fi
                proteins[a.uniprot].add_pathway(pathways[a.pathway])
                pathways[a.pathway].add_protein(proteins[a.uniprot])
            #efor

            print('Adding metabolites')
            metabolites = { p : self.ReactomeNode(p, 'metabolite') for p in pathway_metabolite_map.chemebi }
            for i, a in pathway_metabolite_map.iterrows():
                if a.pathway not in pathways:
                    pathways[a.pathway] = self.ReactomePathway(a.pathway, 'Unknown pathway', a.organism)
                #fi
                metabolites[a.chemebi].add_pathway(pathways[a.pathway])
                pathways[a.pathway].add_metabolite(metabolites[a.chemebi])
            #efor

            print('Adding hierarchy')
            for i, row in self.pathway_hierarchy.iterrows():
                if (row.parent not in pathways) or (row.child not in pathways):
                    continue
                #fi

                parent = pathways[row.parent]
                child  = pathways[row.child]

                parent.add_child(child)
                child.add_parent(parent)
            #efor

            nt = namedtuple('ReactomeIndex', ['pathways', 'proteins', 'metabolites'])

            return nt(pathways, proteins, metabolites)
        #edef
        
        self._obj.register('index', [], lambda f: make_index())
    #edef
    
    ##################################################################
    
    class ReactomePathway(object):
        """
        Internal representation of a reactome pathway.
        
        """
        def __init__(self, identifier, description, organism):
            self._id   = identifier
            self._desc = description
            self._org  = organism
            self._proteins    = []
            self._metabolites = []
            self._children    = []
            self._parents     = []
        #edef
        
        @property
        def id(self):
            """
            Return the identifier of this pathway
            """
            return self._id
        #edef
        
        @property
        def description(self):
            """
            Return a description of this pathway
            """
            return self._desc
        #edef
        
        @property
        def organism(self):
            """
            Return the organism of this pathway
            """
            return self._org
        #edef
        
        @property
        def proteins(self):
            """
            Return the proteins annotated to this object
            """
            return self._proteins
        #edef
        
        @property
        def metabolites(self):
            """
            Return the metabolites annotated to this object
            """
            return self._metabolites
        #edef
        
        @property
        def parents(self):
            """
            Return the parent pathway(s) of this object
            """
            return self._parents
        #edef
        
        @property
        def children(self):
            """
            Return the children pathways of this object
            """
            return self._children
        #edef
        
        @property
        def leaves(self):
            """
            Returns the number of (proteins, metabolites) for all the leaves in the hierarchy from this node onwards.
            """
            nodes = [ self ]
            total_proteins    = 0
            total_metabolites = 0
            for node in nodes:
                if len(node.children) == 0:
                    total_proteins += len(node.proteins)
                    total_metabolites += len(node.metabolites)
                else:
                    nodes.extend(node.children)
                #fi
            #efor
            return (total_proteins, total_metabolites)
        #edef
        
        def add_protein(self, protein):
            """
            Add a protein to this pathway
            parameters:
            -----------
            protein: ReactomeNode object
            """
            if protein not in self._proteins:
                self._proteins.append(protein)
            #fi
        #edef
        
        def add_metabolite(self, metabolite):
            """
            Add a metabolite to this pathway
            parameters:
            -----------
            metabolite: ReactomeNode object
            """
            if metabolite not in self._metabolites:
                self._metabolites.append(metabolite)
            #fi
        #edef
        
        def add_child(self, pathway):
            """
            Add a child pathway to this pathway
            parameters:
            -----------
            pathway: ReactomePathway object
            """
            self._children.append(pathway)
        #edef
        
        def add_parent(self, pathway):
            """
            Add a parent pathway to this pathway
            parameters:
            -----------
            pathway: ReactomePathway object
            """
            self._parents.append(pathway)
        #edef
        
        def __str__(self):
            """
            String representation of this object
            """
            dstr = 'Reactome Pathway\n'
            dstr += ' ID: %s\n' % self.id
            dstr += ' Description: %s\n' % self.description
            dstr += ' Proteins: (%d)\n' % len(self.proteins)
            dstr += ' Metabolites: (%d)\n' % len(self.metabolites)
            dstr += ' Parents: %s' % ', '.join([ p.id for p in self.parents ])
            dstr += ' Children: (%d)\n' % len(self.children)
            prot_leaves, metabolite_leaves = self.leaves
            dstr += ' Leaves (P: %d, M: %d)\n' % (prot_leaves, metabolite_leaves)
            return dstr
        #edef
        
        def __repr__(self):
            """
            String representation of this object
            """
            return str(self)
        #edef
    #eclass
    
    ##################################################################
    
    class ReactomeNode(object):
        """
        A class to handle nodes (proteins/metabolites) in the reactome database
        """
        def __init__(self, identifier, type):
            self._id   = identifier
            self._type = type
            self._pathways = []
        #edef
        
        @property
        def id(self):
            """
            Return the identifier of this object
            """
            return self._id
        #edef
        
        @property
        def type(self):
            """
            Return the type of this object
            """
            return self._type
        #edef
        
        @property
        def pathways(self):
            """
            Return a list of all pathways that this object is annotated to
            """
            return self._pathways
        #edef
        
        def add_pathway(self, pathway):
            """
            Add a pathway to this object
            parameters:
            ----------
            pathway: A ReactomePathway object
            """
            if pathway not in self._pathways:
                self._pathways.append(pathway)
            #fi
        #edef
        
        def __str__(self):
            """
            String representation of this object
            """
            dstr = 'Reactome %s Node\n' % self.type
            dstr += ' ID: %s\n' % self.id
            dstr += ' Pathways: (%d) %s\n' % ( len(self.pathways), ', '.join([p.id for p in self.pathways]))
            return dstr
        #edef
        
        def __repr__(self):
            """
            String representation of this object
            """
            return str(self)
        #edef
    #eclass
            
    
    @property
    def pathways(self):
        """
        Return a list of all pathways identifiers indexed in the reactome database
        """
        return list(self.index.pathways.keys())
    #edef
    
    @property
    def proteins(self):
        """
        Return a list of all protein identifiers indexed in the reactome database
        """
        return list(self.index.proteins.keys())
    #edef
    
    @property
    def metabolites(self):
        """
        Return a list of all metabolites identifiers indexed in the reactome database
        """
        return list(self.index.metabolites.keys())
    #edef
    
    @property
    def organisms(self):
        """
        Return a list of all organisms indexed in this instance of the reactome database
        """
        return self._organisms
    #edef
    
    @property
    def pathway(self):
        """
        Return a dictionary of all pathways indexed in the Reactome database
        
        example:
        
        r = Reactome('Homo sapiens')
        r.protein['R-HSA-8956319']
        """
        return self.index.pathways
    #edef
    
    @property
    def protein(self):
        """
        Return a dictionary of all proteins indexed in the Reactome database
        
        example:
        
        r = Reactome('Homo sapiens')
        r.protein['O43598']
        """
        return self.index.proteins
    #edef
    
    @property
    def metabolite(self):
        """
        Return a dictionary of all metabolites indexed in the Reactome database
        
        example:
        
        r = Reactome('Homo sapiens')
        r.metabolites['warfarin [cytosol]']
        """
        return self.index.metabolites
    #edef
    
    def enrich(self, your_set, background=None, pathway=None, abcd_values=False, method=None, **kwargs):
        """
        Enrich: Check enrichment of Reactome pathways in a given set

        Inputs:
          - your_set: List of Uniprot IDs to test
          - pathway: List of pathways (or single pathway) to test (Defaults to all pathways that your geneIDs are present in)
          - background: List of Entrez Gene IDs to use as background (e.g. set of all expressed genes)
          - abcd_values: Boolean. If True, it will return the actual element values in the contingency table, rather than just counts
          - method: Type of multiple testing correction procedure to use
          - **kwargs: Additional arguments for biu.stats.p_adjust

        Outputs:
         - df : Pandas Data Frame of test results
        """
        your_set = set([ str(ID) for ID in your_set if ID is not None]) & set(self.proteins)

        if pathway is None:
            pathway = set([ p.id for prot in your_set for p in self.protein[prot].pathways if prot in self.proteins ])
        #fi

        if background is None:
            background = self.proteins
        #fi
        background = set(background)

        if isinstance(pathway, str):
            pathway = [ pathway ]
        #fi
        R = []
        for p in pathway:
            pathway_genes = set([prot.id for prot in self.pathway[p].proteins ]) & background
            res = stats.enrichment.setEnrichment(your_set, pathway_genes, background, abcd_values=abcd_values)
            R.append((p, self.pathway[p].id, res.method, res.c2statistic, res.oddsratio, res.pvalue, res.table))
        #efor

        df = pd.DataFrame(R, columns=['pathway', 'name', 'method', 'c2statistic', 'oddsratio', 'p', 'table'])
        if method is not None:
            df['q'] = stats.p_adjust(df.p.values, method, **kwargs)
        #fi

        return df
    #edef
    
    def summary(self, your_set):
        """
        For a set of proteins, count the occurrence of each pathway.

        parameters:
        -----------
        your_set: List[String] A list of uniprot protein identifiers

        Returns:
        Pandas Dataframe of pathway counts
        """

        pathways = [ pway for prot in your_set for pway in self.protein[prot].pathways ]
        pathways = list(ops.lst.freq(pathways).items())
        pathways = pd.DataFrame(pathways, columns=['pathway', 'count'])
        pathways['description'] = pathways.pathway.apply(lambda p: p.description)
        pathways['total_count'] = pathways.pathway.apply(lambda p: len(p.proteins))
        pathways.pathway = pathways.pathway.apply(lambda p: p.id)

        return pathways
    #edef
    
#eclass