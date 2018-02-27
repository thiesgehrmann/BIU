from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings
from .. import utils

###############################################################################

versions = {
  "current" :
  { "vcfCodingURL"     : None,
    "vcfNoncodingURL"  : None,
    "exprURL"          : None,
    "methylURL"        : None
  }
}

def urlFileIndex(version, username, password):

  def sftpCommand(location):
    return "echo -en  'open sftp://sftp-cancer.sanger.ac.uk\\nuser \"%s\" \"%s\"\\ncat \"%s\"' | lftp | zcat" % (username, password, location)
  #edef

  files = { }
  files["vcfCoding"] = (None, 'vcfCoding.vcf.bgz', { "curlCommand" : sftpCommand("cosmic/grch38/cosmic/v84/VCF/CosmicCodingMuts.vcf.gz"),
                                                     "tabix" : True,
                                                     "bgzip" : True,
                                                     "seqField":1,
                                                     "beginField":2, 
                                                     "endField":2 })
  files["vcfCoding_tbi"] = (None, 'vcfCoding.vcf.bgz.tbi', {})
  files["vcfNonCoding"] = (None, 'vcfNonCoding.vcf.bgz', { "curlCommand" : sftpCommand("cosmic/grch38/cosmic/v84/VCF/CosmicNonCodingVariants.vcf.gz"),
                                                           "tabix" : True,
                                                           "bgzip" : True,
                                                           "seqField":1,
                                                           "beginField":2,
                                                           "endField":2 })
  files["vcfNonCoding_tbi"] = (None, 'vcfNonCoding.vcf.bgz.tbi', {})

  return { k : (u, 'cosmic_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for version in versions:
    print(" * %s" % version)
  #efor
#edef

class Cosmic(fm.FileManager):

  def __init__(self, username, password, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version, username, password), objects=[ "vcfCoding", "vcfNonCoding" ], **kwargs)
    self.version = version

    # Define the objects in the clinVar fields
    self.vcfCoding = rm.VCFResourceManager(self, "vcfCoding", "vcfCoding_tbi")
    self.vcfNonCoding = rm.VCFResourceManager(self, "vcfNonCoding", "vcfNonCoding_tbi")

    self.addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  ###############################################################################

  def queryVCF(self, chromosome, start, end):
    return self.vcf.query(chromosome, start, end)
  #edef

  def querySummary(self, chromosome, start, end):
    return self.summary.query(chromosome, start, end, namedtuple=True)
  #edef

#eclass

#eclass
###############################################################################
