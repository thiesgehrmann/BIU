from ..structures import Dataset2
from .. import formats
from .. import utils

np = utils.py.loadExternalModule('numpy')

###############################################################################

class CADD2(Dataset2):
    """
    Load the CADD Dataset.
    
    Example usage
    cadd = biu.db.CADD2()
    cadd.query(12, 500010, 500012)   # Get the CADD scores for a specific region
    cadd.region_thresh(12, 500010, 500012) # Get the 95% of CADD scores in a region.
    """

    versions = { "GRCh37" : {
            "cadd_url"     : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz",
            "cadd_tab_url" : "http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi"
        }
    }

    version = None

    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        """
        Initialize the CADD data structure.
        
        version: Only one version exists at the moment.
        *pargs, **kwargs. See arguments for biu.structures.Dataset2
        
        """
        super(CADD2, self).__init__("CADD/%s" % version, *pargs, **kwargs)
        
        vData = self.versions[version]
        self._obj.add_file('scores.tsv.bgz', utils.Acquire2().curl(vData["cadd_url"]))
        self._obj.add_file('scores.tsv.bgz.tbi', utils.Acquire2().curl(vData["cadd_tab_url"]))
      
        cadd_fields = [ "chrom", "pos", "ref", "alt", "rawscore", "phred" ]
        self._obj.register("scores", [ "scores.tsv.bgz", "scores.tsv.bgz.tbi" ],
                             lambda f: formats.Tabix(f["scores.tsv.bgz"], fieldNames=cadd_fields))
        self.version = version

        self._add_str_func(lambda s: "Version: %s" % self.version)
    #edef

    #############################################################################

    def query(self, chrom, pos, end=None, alt=None):
        """
        Get the CADD scores for a specific region, base, or variant.
        
        parameters:
        -----------
        chrom: String. Which chromosome to look at
        pos:   Integer. What position to look at
        end:   If specified, search from pos until this position
        alt:   If specified, return only CADD scores for this alternative base.
        
        returns:
            if end is None:
                
            if end is None and alt is not None:
                Returns a float
            Otherwise:
                returns a dictionary of alternative bases and cadd scores.
        """
        qres = self.scores.query(chrom, pos, (pos if (end is None) else end), namedtuple=True)
        if (alt is None) and (end is None):
            resPhred = {}
            for res in qres:
                resPhred[res.alt] = float(res.phred)
            #efor
            return resPhred
        #fi
        if (alt is None) and (end is not None):
            resPhred = {}
            for res in qres:
                resPhred[chrom, (int(res.pos),res.alt)] = float(res.phred)
            #efor
            return resPhred
        #fi
        if (alt is not None) and (end is not None):
            resPhred = {}
            for res in [r for r in qres if r.alt == alt]:
                resPhred[(chrom, int(res.pos))] = float(res.phred)
            #efor
            return resPhred
        #fi
        if (end is None) and (alt is not None):
            relRes = [ r for r in qres if r.alt == alt ]
            if len(relRes) != 1:
                return None
            else:
                return float(relRes[0].phred)
            #fi
        #fi
    #edef

    def query_regions(self, regions):
        """
        Get the CADD scores for multiple regions.
        
        parameters:
        -----------
        regions: List of 3-tuples of (chromosome, start, end), or 2-tuples of (chromosome, pos)
        
        returns:
        a dict of cadd scores.
        """
        resPhred = {}
        for reg in regions:
            c = None
            s = None
            e = None
            if len(reg) == 3:
                c,s,e = reg
            elif len(reg) == 2:
                c,s = reg
            else:
                raise ValueError("Incorrect region specified")
            #fi
            
            qres = self.scores.query(c, s, e, namedtuple=True) 
            for res in qres:
                resPhred[(c, int(res.pos),res.alt)] = float(res.phred)
            #efor
        #efor
        return resPhred
    #edef

    def region_thresh(self, chrom, start, end, percentile=95):
        """
        Determine the 95% threshold of CADD scores for a given region
        
        parameters:
        -----------
        chrom: String. Which chromosome to look at
        pos:   Integer. What position to look at
        end:   If specified, search from pos until this position
        percentile: what percentile to use.

        returns:
        A float
        
        """
        return self.regions_thresh( [(chrom, start, end)], percentile )
    #edef

    def regions_thresh(self, regions, percentile=95):
        """
        parameters:
        -----------
        regions: List of 3-tuples of (chromosome, start, end), or 2-tuples of (chromosome, pos)
        percentile: what percentile to use.
        
        returns:
        A float
        """
        return np.percentile(np.array(list(self.query_regions(regions).values())), percentile)
    #edef

#eclass