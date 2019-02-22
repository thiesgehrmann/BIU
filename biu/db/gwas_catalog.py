from .. import ops
from .. import utils
from .. import settings
from ..structures import Dataset2

pd = utils.py.loadExternalModule("pandas")
plt = utils.py.loadExternalModule("matplotlib", "pylab")

###############################################################################

class GWAS_Catalog(Dataset2):
    """
    An interface to the GWAS Catalog https://www.ebi.ac.uk/gwas
    
    Available in two versions:
     * 1.0.2
     * 1.0
    
    Provides several objects:
    
     * assoc : A table of associations listed in the GWAS catalog
     * studies: A table of the studies
     * ancestry: A table of the ancestry information in the studies
     
    Also provides several functions:
    
     * query: Return only associations that pass certain filters
     * plot:  Plot the associations
    """

    version = None

    versions = {
        "1.0.2" : {
          "assoc" : "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
          "studies" : "https://www.ebi.ac.uk/gwas/api/search/downloads/studies",
          "ancestry" : "https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry"
        },
        "1.0" : {
          "assoc" : "https://www.ebi.ac.uk/gwas/api/search/downloads/full",
          "studies" : "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative",
          "ancestry" : "https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry"
        }
    }

    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        """
        Initialize the GWAS_Catalog object.
        
        parameters:
        -----------
        version: String. [ "1.0.2", "1.0"]
        
        for *pargs and **kwargs see the docstring for Dataset2
        """
        super(GWAS_Catalog, self).__init__("GWAS_Catalog/%s" % version, *pargs, **kwargs)
        self.version = version
        
        self._add_str_func( lambda s: "Version : %s" % s.version )
        
        self._obj.add_file('assoc.tsv',    utils.Acquire2().curl(self.versions[version]["assoc"]))
        self._obj.add_file('studies.tsv',  utils.Acquire2().curl(self.versions[version]["studies"]))
        self._obj.add_file('ancestry.tsv', utils.Acquire2().curl(self.versions[version]["ancestry"]))
        
        def loadAssoc(d):
            assoc = pd.read_csv(d["assoc.tsv"], sep='\t')

            assoc = assoc[~assoc.CHR_POS.apply(lambda x: 'x' in str(x))]

            grouped   = assoc.CHR_POS.apply(lambda x: ';' in str(x))
            ungrouped = ops.dataframe.flat(assoc[grouped], sep=';', ignore_unequal_rows=True,
                                           fields=[ 'CHR_ID', 'CHR_POS', 'MAPPED_GENE', 'SNPS',
                                                    'CONTEXT', 'STRONGEST SNP-RISK ALLELE' ])

            total = pd.concat([assoc[~grouped], ungrouped])

            field_types = [ ('P-VALUE', float, 1),
                            ('CHR_POS', int, 0) ]

            for (field, field_type, nan_value) in field_types:
                total[field] = total[field].fillna(nan_value).astype(field_type)
            #efor

            return total
        #edef
        
        self._obj.register( "assoc", ["assoc.tsv"], loadAssoc)
        
        self._obj.register("studies",  ["studies.tsv"],
                           lambda d: pd.read_csv(d["studies.tsv"], sep='\t'))
        
        self._obj.register("ancestry", ["ancestry.tsv"],
                           lambda d: pd.read_csv(d["ancestry.tsv"], sep='\t', index_col=False))

    #edef
    
    def plot(self, chrom, start=None, stop=None, pvalue=None, trait=None, ax=None, **scatter_args):
        """
        Plot the GWAS hits in a specific region.

        Parameters:
        -----------
        chrom:  String.  The chromosome ID to use
        start:  Integer. The chromosomal position to start...
        stop:   Integer. and end at.
        pvalue: Float.   Only tests with pvalues lower than this value
        trait:  String|List of strings. Only tests in these given traits.
        ax:     Matplotlib axis. The axis to plot onto
        **scatter_args: dict. Additional arguments to scatter for plotting.

        Output:
        -------

        Tuple of: (Matplotlib figure, Matplotlib axis plotted to)
        """

        if ax is None:
            fig, axes = utils.figure.subplots(ncols=1, nrows=1, figsize=(10, 5))
            ax = axes[0]
        #fi

        rel = self.query(chrom, start, stop, pvalue, trait)

        ax.scatter(rel.CHR_POS, -np.log10(rel['P-VALUE']), **scatter_args)
        ax.set_ylabel('-log10 pvalue')
        ax.set_xlabel('Genomic Position (Chromosome %s)' % str(chrom))

        return ax.get_figure(), ax
    #edef
    
    def query(self, chrom, start=None, stop=None, pvalue=None, trait=None):
        """
        Query the associations in a given region.
        Parameters:
        -----------
        chrom:  String.  The chromosome ID to use
        start:  Integer. The chromosomal position to start...
        stop:   Integer. and end at.
        pvalue: Float.   Only return pvalues lower than this value
        trait:  String|List of strings. Only tests in these given traits.

        
        
        Output:
        --------
        Pandas DataFrame
        """
        
        if trait is not None:
            if isinstance(trait, str):
                trait = [trait]
            #fi
        #fi
        
        rel = self.assoc[self.assoc.CHR_ID == str(chrom)]
        rel = rel if start is None else rel[rel.CHR_POS >= start]
        rel = rel if stop  is None else rel[rel.CHR_POS <= stop]
        rel = rel if pvalue is None else rel[rel['P-VALUE'] <= pvalue]
        rel = rel if trait is None else rel[rel['DISEASE/TRAIT'].isin(trait)]
        
        return rel
    #edef
        
#eclass