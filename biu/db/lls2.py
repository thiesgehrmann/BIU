from ..structures import Dataset2
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule("pandas")

###############################################################################

class LLS2(Dataset2):
    """
    The Leiden Longevity study object.
    Provides functionality to query the genetic variants in the study.
    
    objects:
    vcf_*, where * is a chromosome identifier, 1-22,M,X,Y
    phenotypes: A description of phenotypes for each individual in the study
    ids: An ID mapping for individuals in the study
    percentiles: Data on the percentiles of ageing for the individuals in the study
    
    Example usage:
    ---------------
    
    lls = biu.db.LLS2()
    v = lls.filter(15, 5001234, 5002345)
    
    """

    versions = { "current":
        { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "M", "X", "Y" ] }
    }
    version = None
    vcf = None

    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        """
        Initialize the objet.
        
        See parameters for biu.structures.Dataset2 for details about *pargs and **kwargs
        """
        super(LLS2, self).__init__("LLS/%s" % version, *pargs, **kwargs)
        
        self.version = version
        
        for chrom in self.versions[self.version]["chrs"]:
            vcf = "vcf_%s.vcf" % chrom
            tbi = "vcf_%s.vcf.tbi" % chrom
            vcf_file = utils.Acquire2("/exports/molepi/LLSSEQ/tbx/merged.chr%s.vcf.bgz" % chrom)
            tbi_file = utils.Acquire2("/exports/molepi/LLSSEQ/tbx/merged.chr%s.vcf.bgz.tbi" % chrom)
            self._obj.add_file(vcf, vcf_file, finalize=False)
            self._obj.add_file(tbi, tbi_file, finalize=False)

            self._obj.register("vcf_%s" % chrom, [vcf, tbi], lambda f, vcf=vcf: formats.VCF2(f[vcf], tabix=True))
        #efor
        
        self._obj.add_file('percentiles', utils.Acquire2("/exports/molepi/tgehrmann/data/LLS_data/lls_percentiles.tsv"), finalize=False)
        self._obj.add_file('ids', utils.Acquire2('/exports/molepi/tgehrmann/data/LLS_data/lls_ids.csv'), finalize=False)
        self._obj.add_file('phen', utils.Acquire2("/exports/molepi/tgehrmann/data/LLS_data/lls_phenotypes218.txt"), finalize=False)

        self._obj.register('percentiles', ['percentiles'], lambda f: pd.read_csv(f[ "percentiles" ], sep='\t'))
        self._obj.register('ids',         ['ids'],         lambda f: pd.read_csv(f[ "ids"]))
        
        def load_phen(f):
            dat = pd.read_csv(f["phen"], delimiter=' ')
            dat['famnr'] = dat.LLnr.apply(lambda x: str(int(x.split('.')[0])))
            return dat
        #edef
        
        self._obj.register('phenotypes', ['phen'], load_phen)
    
        self._add_str_func(lambda s: "Version: %s" % self.version)
    #edef


    def filter(self, chrom, start, end, *pargs, **kwargs):
        """
        Perform a filter for several regions
        
        parameters:
        -----------
        chrom: String|int. Chromosome of interest
        start: int. Start of region of interest
        end: int. End of region of interest
        *pargs, **kwargs: See additional arguments for VCF2.filter
        
        Returns: VCF2 object
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi
        
        return self._obj[oname].filter(chrom, start, end, *pargs, **kwargs)
    #edef
    
    def filter_regions(self, regions, chrom=None, start=None, end=None, *pargs, **kwargs):
        """
        Perform a filter for several regions
        
        parameters:
        -----------
        regions: A list of 3-tuples (chrom, start, end) for each region of interest
        chrom, start, end: Ignored
        *pargs, **kwargs: See additional arguments for VCF2.filter
        
        Returns: VCF2 object
        """
        rets = [ self.filter(c, s, e, *pargs, **kwargs) for (c,s,e) in regions ]
        return rets[0].merge(rets)
    #edef


    def get_var(self, chrom, *pargs, **kwargs):
        """
        Get the variant record for a specific variant
        
        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?
        
        returns:
        --------
        A cyvcf2.VCF.Variant object if the variant exists. Otherwise None
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi

        return self._obj[oname].get_var(chrom, *pargs, **kwargs)
    #edef

    def who_has(self, chrom, *pargs, **kwargs):
        """
        Determine who has a specific variant

        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?

        returns:
        --------
        List of sample IDs for who has the variant.
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi

        return self._obj[oname].who_has(chrom, *pargs, **kwargs)
    #edef


    def cgid2llnr(self, cgID):
        """
        Convert a cgID (sequencing identifier) to LLnr (study number)
        
        parameters:
        -----------
        cgID: String. The sequencing identifier
        
        returns: LLnr string.
        """
        cgID    = cgID.replace('_240_37-ASM', '')
        possible = self.phenotypes[self.phenotypes.cgID == cgID].LLnr.values 
        if len(possible) > 0:
            return possible[0]
        #fi
        return None
    #edef

    def llnr2cgid(self, llnr):
        """
        Convert a LLNR (study number) to cgID  (sequencing identifier)
        
        parameters:
        -----------
        llnr: String. The LLnr
        
        returns: cgID
        """
        possible = self.phenotypes[self.phenotypes.LLnr == llnr].cgID.values
        if len(possible) > 0:
            return possible[0]
        #fi
        return None
    #edef

    def pass_fams(self, pThresh=.90):
        """
        Determine which families pass the novel selection threshold.
        I.e. at least two siblings > pThresh %, and at least one parent > pThresh %

        Inputs:
          pThresh: Float. Percentile threshold to satisfy
        Output:
          Set of family IDs (Integers)
        """
        p = self.percentiles
        pF0 = p[(p.percentile >= pThresh) & (p.generation == 1)]
        pF1 = p[(p.percentile >= pThresh) & (p.generation == 2)]

        f0Pass = pF0.famnr.unique()
        f1Pass = pF1[ pF1.groupby("famnr").transform(len).llnr >= 2 ].famnr.unique()

        famPass = set(f0Pass) & set(f1Pass)
        return famPass
    #edef

    def pass_indiv(self, pThresh=.90, relFams=None, sequenced=False):
        """
        Determine, in the families that pass the novel selection threshold, which individuals satisify the criteria.
        Inputs:
          pThresh: Float. Percenfile threshold to satisfy
          relFams: Set of Integers. Family IDs to consider (default is None, then output of pass_fams is used)
          sequenced: Boolean. Only return IDs of individuals that have been sequenced. (Default False)
        Output:
          Set of individual IDs.
        """
        # Select all individuals from that family

        if relFams is None:
            relFams = self.pass_fams(pThresh)
        #fi

        pPer = self.percentiles
        relIndivs = pPer[pPer.famnr.apply(lambda f: f in relFams) & (pPer.percentile >= pThresh)]
        relIndivs = relIndivs[["famnr", "llnr", "percentile", "motherID", "fatherID", "gender"]]

        relIndiv = set(relIndivs.llnr.values)
        if sequenced:
            relIndiv = relIndiv & set(self.phenotypes.LLnr.values)
        #fi

        return relIndiv
    #edef


#eclass
