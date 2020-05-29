from ..structures import Dataset2
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule("pandas")
import numpy as np

from functools import reduce

###############################################################################

class VCF_gnomad(formats.VCF2):
    def genotype_matrix(self):
        raise NotImplementedError("This function is not available for the GNOMAD dataset.")
    #edef
#eclass

class GnomadBiuVariant(formats.BiuVariant):
     def summary(self, altpos=None, sub=[None,'AFR','AMR','ASJ','EAS','FIN','NFE','OTH','SAS']):
        """
        Make a summary of the variants at this variant
        ASSUMES A MONOPLOID/DIPLOID ORGANISM!
        parameters:
        self: BiuVariant Variant object
        altpos: Integer, List[Integer]
            Which alternative allele indexes to consider
        sub: List[String]
            Which subpopulations to consider:
            None : World
            'AFR'
            'AMR'
            'ASJ'
            'EAS'
            'FIN'
            'NFE'
            'OTH'
            'SAS'
            
        Returns: DataFrame
        """
        
        counts  = []
        idxs = []
        
        if altpos is None:
            altpos = list(range(len(self.ALT)))
        elif isinstance(altpos, int):
            altpos = [altpos]
        #fi
        
        def mk_summary(var, alt_pos, sub):
            ID  = self.make_identifier(alt_pos)
            dat = []
            for sub_name in sub:
                if sub_name is None:
                    dat = dat + all_summary(var, alt_pos)
                else:
                    dat = dat + sub_summary(var, alt_pos, sub_name)
                #fi
            #efor
                    
            return [ID] + dat
        #edef
            
        
        def all_summary(var, altp):
            gcIndexes = formats.VCF2.genotype_info_field_indexes(altp+1)
            
            gcmale    = [ var.INFO.get("GC_Male")[i] for i in gcIndexes ] if (var.INFO.get("GC_Male") is not None) else [0, 0, 0]
            gcfemale  = [ var.INFO.get("GC_Female")[i] for i in gcIndexes ] if (var.INFO.get("GC_Female") is not None) else [0, 0, 0]
            
            rr = 0
            r  = 0
            ra = 0
            a  = 0
            aa = 0
            u  = 0
            
            chrom = var.CHROM.lower()     
            total = var.INFO.get("AN_Male") + var.INFO.get("AN_Female")
            if chrom == 'y':
                a = var.INFO.get("AC")[altp]
                r = var.INFO.get("AN") - a
            elif chrom == 'x':
                rr = gcfemale[0]
                ra = gcfemale[1]
                aa = gcfemale[2]
                r = gcmale[0]
                a = gcmale[1]
            else:
                rr = gcmale[0] + gcfemale[0]
                ra = gcmale[1] + gcfemale[1]
                aa = gcmale[2] + gcfemale[2]
            #fi
            
            u = total - (2*rr + r + 2*ra + a + 2*aa)
            af = var.INFO.get("AF")[altp] if len(var.ALT) > 1 else var.INFO.get("AF")

            return [ rr, r, ra, a, aa, u, af ]
        #edef

        def sub_summary(var, altp, sub):
            # See http://gnomad.broadinstitute.org/faq for the different subpopulations

            gc = "GC_%s" % sub
            gcIndexes = formats.VCF2.genotype_info_field_indexes(altp+1)
            gc        = [ var.INFO[gc][i] for i in gcIndexes ] if (var.INFO.get(gc) is not None) else [0, 0, 0]
            
            total = var.INFO.get('AN_%s' % sub)
            
            rr = gc[0]
            r  = 0
            ra = gc[1]
            a  = 0
            aa = gc[2]
            u  = total - 2*(rr+ra+aa)
            
            af = var.INFO.get("AF_%s" % sub)[altp] if len(var.ALT) > 1 else var.INFO.get("AF_%s" % sub)
            
            return [rr, r, ra, a, aa, u, af]
        #edef
        
        S = [ mk_summary(self, alt_pos, sub) for alt_pos in range(len(self.ALT)) ]
        
        sub_names = [ 'ALL' if s is None else s for s in sub ]
        cols = ["index"] + [ "%s_%s" % (sub, allele)
                         for sub in sub_names
                             for allele in  ["RR", "R", "RA", "A", "AA", "O", "AF"]
                        ]
        S = pd.DataFrame(S, columns=cols).set_index("index")
        return S
    #edef
#eclass
        

class Gnomad2(Dataset2):
    """
    An interface to the GNOMAD dataset
    """

    versions = { "GRCh37_2.0.2" : {
        "vcf"             : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz",
        "tbi"             : "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi",
        "cov_url_proto"   : "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr%s.coverage.txt.gz",
        "chr"             : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" ],
        "cov_seq_field"   : 1,
        "cov_begin_field" : 2,
        "cov_end_field"   : 2
        }
    }

    version = None

    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        super(Gnomad2, self).__init__("Gnomad/%s" % version, *pargs, **kwargs)
        self.version = version
        
        v_data = self.versions[version]
    
        self._obj.add_file("var.vcf.bgz", utils.Acquire2().curl(v_data['vcf']))
        self._obj.add_file("var.vcf.bgz.tbi", utils.Acquire2().curl(v_data['tbi']))
        self._obj.register("vcf", ["var.vcf.bgz", "var.vcf.bgz.tbi"],
                           lambda f: VCF_gnomad(f["var.vcf.bgz"], tabix=True,
                                               internal_variant_representation=GnomadBiuVariant))
        
        cov_entry_fields = [ "chrom", "pos", "mean", "median", "q1", "q5", "q10", "q15", "q20", "q25", "q30", "q50", "q100" ]
    
        for chrom in self.versions[version]["chr"]:
            cov     = utils.Acquire2().curl(v_data["cov_url_proto"] % chrom).gunzip().bgzip()
            cov_tbi = cov.tabix(seq=v_data['cov_seq_field'], start=v_data['cov_begin_field'], end=v_data['cov_end_field'])
            
            cov_file     = "cov_%s.bgz" % chrom
            cov_tbi_file = "cov_%s.bgz.tbi" % chrom
            
            self._obj.add_file(cov_file, cov)
            self._obj.add_file(cov_tbi_file, cov_tbi)
            
            self._obj.register("cov_%s" % chrom, [cov_file, cov_tbi_file],
                               lambda f, n=cov_file: formats.Tabix(f[n], fieldNames=cov_entry_fields))
        #efor
    
        self._add_str_func(lambda s: "Version: %s" % self.version)
    #edef

  #############################################################################

    def filter(self, *args, **kwargs):
        """
        Perform a filter for several regions

        parameters:
        -----------
        *pargs, **kwargs: See arguments for VCF2.filter

        Returns: VCF2 object
        """
        return self.vcf.filter(*args, **kwargs)
    #edef

    def filter_regions(self, regions, chrom=None, start=None, end=None, *args, **kwargs):
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
        return reduce(lambda a,b: a+b, rets)
    #edef

    def filter_cov(self, chrom, start, end, **kwargs):
        """
        Get the coverage for specific regions in the GNOMAD dataset
        parameters:
        -----------
        chromosome: String. Chromosome
        start:      int. Start of region of interest
        end:        int. End of region of interest
        **kwargs: See parameters for biu.formats.Tabix.query
        
        returns:
        --------
        A list of entries in the coverage file 
        """
        chrom = str(chrom)
        oname = "cov_%s" % chrom
        return self._obj.get(oname).query(chrom, start, end, **kwargs)
    #edef

    def filter_cov_regions(self, regions, **kwargs):
        """
        Get the coverage for specific regions in the GNOMAD dataset
        parameters:
        -----------
        regions: list of 3-tuples of (chrom, start, end).
        chromosome: String. Chromosome
        start:      int. Start of region of interest
        end:        int. End of region of interest
        **kwargs: See parameters for biu.formats.Tabix.query
        
        returns:
        --------
        A list of entries in the coverage file 
        """
        R = []
        for (c,s,e) in regions:
            for r in self.query_cov(c, s, e, **kwargs):
                R.append(r)
            #efor
        #efor
        return R
    #edef 

    def who_has(self, *pargs, **kwargs):
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
        return self.vcf.who_has(*pargs, **kwargs)
    #edef

    def get_var(self, *pargs, **kwargs):
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
        return self.vcf.get_var(*pargs, **kwargs)
    #edef

#eclass
