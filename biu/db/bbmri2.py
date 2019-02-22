from ..structures import Dataset2
from .. import formats
from .. import utils

###############################################################################

class BBMRI2(Dataset2):
    """
    The BBMRI complete genomics sequencing set.
    Provides functionality to query the genetic variants in the study.
    
    objects:
    vcf_*, where * is a chromosome identifier, 1-22,M,X,Y
    phenotypes: A description of phenotypes for each individual in the study
    ids: An ID mapping for individuals in the study
    percentiles: Data on the percentiles of ageing for the individuals in the study
    
    Example usage:
    ---------------
    
    lls = biu.db.BBMRI2()
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
        super(BBMRI2, self).__init__("BBMRI/%s" % version, *pargs, **kwargs)
        
        self.version = version
        
        for chrom in self.versions[self.version]["chrs"]:
            vcf = "vcf_%s.vcf" % chrom
            tbi = "vcf_%s.vcf.tbi" % chrom
            vcf_file = utils.Acquire2("/exports/molepi/BBMRISEQ/tbx/merged.bbmri.chr%s.vcf.bgz" % chrom)
            tbi_file = utils.Acquire2("/exports/molepi/BBMRISEQ/tbx/merged.bbmri.chr%s.vcf.bgz.tbi" % chrom)
            self._obj.add_file(vcf, vcf_file, finalize=False)
            self._obj.add_file(tbi, tbi_file, finalize=False)

            self._obj.register("vcf_%s" % chrom, [vcf, tbi], lambda f, vcf=vcf: formats.VCF2(f[vcf], tabix=True))
        #efor
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

#eclass
