from ..structures import Dataset2
from .. import formats
from .. import utils
from .. import ops

import json

class Cosmic(Dataset2):
    """
    Access the COSMIC database
    
    In versions: grch37 or grch38
    
    objects:
    --------
    coding: coding variants
    noncoding: noncoding variants
    
    functions:
    ----------
    query: Perform a query
    queryRegions: Perform multiple queries in one go.
    """
    
    versions = {
        "grch37" : {
            "coding"    : "https://cancer.sanger.ac.uk/cosmic/file_download?data=GRCh37%2Fcosmic%2Fv87%2FVCF%2FCosmicCodingMuts.vcf.gz",
            "noncoding" : "https://cancer.sanger.ac.uk/cosmic/file_download?data=GRCh37%2Fcosmic%2Fv87%2FVCF%2FCosmicNonCodingVariants.vcf.gz"
        },
        "grch38" : {
            "coding"    : "https://cancer.sanger.ac.uk/cosmic/file_download?data=GRCh38%2Fcosmic%2Fv87%2FVCF%2FCosmicCodingMuts.vcf.gz",
            "noncoding" : "https://cancer.sanger.ac.uk/cosmic/file_download?data=GRCh38%2Fcosmic%2Fv87%2FVCF%2FCosmicNonCodingVariants.vcf.gz"
        }
    }
    
    class CosmicAcquire(utils.Acquire2):
        """
        Additional Acquire functionality to cover the crazy COSMIC authentication steps.
        See https://cancer.sanger.ac.uk/cosmic/help/file_download
        """
        def cosmic(self, username, password, url):            
            curl_hash   = ops.lst.hash([ url, username, password ])
            output_file = self.AcquireFile(dirname=None, basename=curl_hash)
            
            def _cosmic_curl(url, username, password, output_file):
                auth = utils.exe.getCommandOutput("echo '%s:%s' | base64" % (username, password), shell=True).decode().strip()
                request_json = utils.exe.getCommandOutput('curl -H "Authorization: Basic %s" "%s"' % (auth, url))
                
                request_json = json.loads(request_json)
                
                request_url = request_json["url"]
                p = utils.exe.runCommand("curl '%s' > '%s'" % (request_url, output_file), shell=True, verbose=True)
                return p
            #edef
            
            step = self.AcquireStep("Cosmic(%s)" % url, [], output_file,
                                    lambda i,o: _cosmic_curl(url, username, password, o))
            return self.add_step(step)
        #edef
    #eclass
    
    def __init__(self, version=list(versions.keys())[0],
                 username="t.gehrmann@lumc.nl", password="Cosmic_password1",
                 *pargs, **kwargs):
        """
        Initialize the cosmic dataset.
        
        parameters:
        -----------
        version:  [grch37|grch38]
        username: Your COSMIC username
        password: Your COSMIC password
        *pargs, **kwargs: See parameters to Dataset2.
        """
        super(Cosmic, self).__init__("Cosmic/%s" % version, *pargs, **kwargs)
        
        coding    = self.CosmicAcquire().cosmic(username, password, self.versions[version]["coding"]).gunzip().bgzip()
        noncoding = self.CosmicAcquire().cosmic(username, password, self.versions[version]["noncoding"]).gunzip().bgzip()

        coding_tbi = coding.tabix(seq=1, start=2, end=2)
        noncoding_tbi = noncoding.tabix(seq=1, start=2, end=2)
        
        self._obj.add_file("coding.vcf.bgz", coding)
        self._obj.add_file("coding.vcf.bgz.tbi", coding_tbi)
        
        self._obj.add_file("noncoding.vcf.bgz", noncoding)
        self._obj.add_file("noncoding.vcf.bgz.tbi", noncoding_tbi)
        
        self._obj.register("coding", ["coding.vcf.bgz", "coding.vcf.bgz.tbi"],
                           lambda i: formats.VCF(i["coding.vcf.bgz"], tabix=True))
        self._obj.register("noncoding", ["noncoding.vcf.bgz", "noncoding.vcf.bgz.tbi"],
                           lambda i: formats.VCF(i["noncoding.vcf.bgz"], tabix=True))
    #edef
    
    def query(self, chrom, start, stop, coding=True, noncoding=True, *pargs, **kwargs):
        """
        Perform a query on a single location.
        
        parameters:
        -----------
        chrom: The chromosome of interest
        start: Where to start looking
        stop: Where to stop looking
        coding|noncoding: Boolean. Specify whether or not you want to have coding/noncoding results.
        *pargs, **kwargs: See options for biu.formats.VCF.query for more details.
        
        Returns a biu.formats.VCF object
        """
        return self.queryRegions( [ tuple(chrom, start, stop) ], *pargs, **kwargs)
    #edef
    
    def queryRegions(self, regions, coding=True, noncoding=True, extract=None, *pargs, **kwargs):
        """
        Perform a query on a multiple locations.
        
        parameters:
        -----------
        regions: List of tuples of (chrom, start, stop)
        coding|noncoding: Boolean. Specify whether or not you want to have coding/noncoding results.
        *pargs, **kwargs: See options for biu.formats.VCF.query for more details.
        
        Returns a biu.formats.VCF object
        """
        cResults  = self.coding.queryRegions(*args, extract='raw', *pargs, **kwargs) if coding else []
        ncResults = self.noncoding.queryRegions(*args, extract='raw', *pargs, **kwargs) if noncoding else []
        return formats.VCF(cResults + ncResults, self.coding.template)
    #edef
#eclass
###############################################################################
