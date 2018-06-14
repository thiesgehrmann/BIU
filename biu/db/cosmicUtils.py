from ..structures import Dataset
from ..config import settings as settings
from .. import formats
from .. import utils
import pandas as pd

from ..formats import VCF


###############################################################################

class Cosmic(Dataset):

  versions = {
    "grch37_85" :
    { "vcfCodingURL"     : "cosmic/grch37/cosmic/v85/VCF/CosmicCodingMuts.vcf.gz",
      "vcfNoncodingURL"  : "cosmic/grch37/cosmic/v85/VCF/CosmicNonCodingVariants.vcf.gz",
    },
    "grch37_84" :
    { "vcfCodingURL"     : "cosmic/grch37/cosmic/v84/VCF/CosmicCodingMuts.vcf.gz",
      "vcfNoncodingURL"  : "cosmic/grch37/cosmic/v84/VCF/CosmicNonCodingVariants.vcf.gz",
    },
    "grch38_84" :
    { "vcfCodingURL"     : "cosmic/grch38/cosmic/v84/VCF/CosmicCodingMuts.vcf.gz",
      "vcfNoncodingURL"  : "cosmic/grch38/cosmic/v84/VCF/CosmicNonCodingVariants.vcf.gz",
    },
    "grch38_85" :
    { "vcfCodingURL"     : "cosmic/grch38/cosmic/v85/VCF/CosmicCodingMuts.vcf.gz",
      "vcfNoncodingURL"  : "cosmic/grch38/cosmic/v85/VCF/CosmicNonCodingVariants.vcf.gz",
    }
  }

  def __init__(self, username="t.gehrmann@lumc.nl", password="Cosmic_password1", version=list(versions.keys())[0], where=None, **kwargs):
    fileIndex = self.__genFileIndex(version, username, password, where)
    Dataset.__init__(self, fileIndex, **kwargs)
    self.version = version

    self._registerObject("vcf_c", formats.VCF, [ 'vcf_coding', 'vcf_coding_tbi' ], fileIndex['vcf_coding'].path, tabix=True)
    self._registerObject("vcf_nc", formats.VCF, [ 'vcf_non_coding', 'vcf_non_coding_tbi' ], fileIndex['vcf_non_coding'].path, tabix=True)

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, username, password, where=None):
     finalPath = '%s/cosmic_%s' % ( (settings.getWhere() if where is None else where), version)
     files = {}
     vcfCoding = utils.Acquire(where=where).lftp("sftp://sftp-cancer.sanger.ac.uk",
                                                 self.versions[version]["vcfCodingURL"],
                                                 username=username, password=password).gunzip().bgzip()
     files['vcf_coding']     = vcfCoding.finalize('%s/vcfCoding.vcf.bgz' % finalPath)
     files['vcf_coding_tbi'] = vcfCoding.tabix(seq=1, start=2, end=2).finalize('%s/vcfCoding.vcf.bgz.tbi' % finalPath)

     vcfNonCoding = utils.Acquire(where=where).lftp("sftp://sftp-cancer.sanger.ac.uk",
                                                    self.versions[version]["vcfNoncodingURL"],
                                                    username=username, password=password).gunzip().bgzip()
     files['vcf_non_coding']     = vcfNonCoding.finalize('%s/vcfNonCoding.vcf.bgz' % finalPath)
     files['vcf_non_coding_tbi'] = vcfNonCoding.tabix(seq=1, start=2, end=2).finalize('%s/vcfNonCoding.vcf.bgz.tbi' % finalPath)

     return files
  #edef

  ###############################################################################

  def queryVCF(self, *args, **kwargs):
    return self.query(*args, **kwargs)
  #edef

  def query(self, *args, **kwargs):
    return self.queryRegions( [ tuple(args) ], **kwargs)
  #edef

  def queryRegions(self, *args, coding=True, noncoding=True, extract=None, **kwargs):
    cResults  = self._getObject("vcf_c").queryRegions(*args, extract='raw', **kwargs) if coding else []
    ncResults = self._getObject("vcf_nc").queryRegions(*args, extract='raw', **kwargs) if noncoding else []
    return VCF(cResults + ncResults, self._getObject("vcf_c").template)
  #edef

#eclass

#eclass
###############################################################################
