
import json

from ..structures import fileManager as fm
from ..structures import resourceManager as rm
from ..config import settings as settings

from .. import utils
from .. import formats

versions = {
  "human_9606_b150_GRCh37p13" : {
    "vcf" : "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz",
    "tbi" : "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi",
    "assemblyid" : "GRCh37.p13"
  },
  "human_9606_b150_GRCh38p7" : {
    "vcf" : "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz",
    "tbi" : "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz.tbi",
    "assemblyid" : "GRCh38.p7"
  }
}

def urlFileIndex(version):
  files = {}
  files["vcf"] = (versions[version]["vcf"], "all_snps.vcf.bgz", {})
  files["tbi"] = (versions[version]["tbi"], "all_snps.vcf.bgz.tbi", {})
  files["info"] = (None, "info_dict.sqldict.sqlite", {})
  return { k : (u, 'dbSNP_%s/%s' % (version, l), o) for (k, (u, l, o)) in files.items() }
#edef

def listVersions():
  print("Available versions:")
  for genome in versions:
    print(" * %s" % genome)
  #efor
#edef

###############################################################################

class DBSNP(fm.FileManager):

  version = None
  where    = None
  fileIndex = None

  def __init__(self, version=list(versions.keys())[0], **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(version), objects=["vcf", "info"], **kwargs)
    self.version = version

    self.addStrFunction( lambda s: "Version : %s" % s.version )
    self.vcf = rm.VCFResourceManager(self, "vcf", "tbi")
    self.info = rm.SQLDictResourceManager(self, "info")

  #edef

  ###############################################################################

  def query(self, *pargs, **kwargs):
    return self.vcf.query(*pargs, **kwargs)
  #edef

  def queryRegions(self, *pargs, **kwargs):
    return self.vcf.queryRegions(*pargs, **kwargs)
  #edef

  def __getitem__(self, c):
    #Figure out a good way to get info about the rsID from the "pos" file.
    # Tabix doesn't seem to work...
    res = None
    if isinstance(c, int):
      res = self.idLookup(c)
    elif isinstance(c, str):
      c = c.lower()
      if (c[:2] == 'rs') or (c[:2] == 'ss'):
        try:
          vartype = c[:2]
          c = int(c[2:])
        except:
          utils.error("Not a valid dbsnp id: '%s'" % c)
          return None
        #etry
      res = self.idLookup(c)
    else:
      return None
    #fi

    return res
  #edef

  def __getAssemblyPosition(self, ID):
    utils.dbm("Querying for rs%d via REST." % ID)
    assemblyid = versions[self.version]["assemblyid"]

    def seqPosLookup(restResult):
      for asm in restResult['primary_snapshot_data']['placements_with_allele']:
        for seq in asm['placement_annot']['seq_id_traits_by_assembly']:
          if seq['assembly_name'] == assemblyid:
            possiblePositions = [ (allele['allele']['spdi']['seq_id'], allele['allele']['spdi']['position']) for allele in asm['alleles'] ]
            return possiblePositions[0]
          #fi
        #efor
      #efor
      return None
    #edef

    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/%d' % ID
    dat = str(utils.getCommandOutput('curl -X GET --header "Accept: application/json" "%s"' % url)[0].decode('UTF-8'))
    dat = json.loads(dat)

    if 'primary_snapshot_data' not in dat:
      newkey = int(dat['merged_snapshot_data']['merged_into'][0])
      utils.dbm("Redirecting %d -> %d" % (ID, newkey))
      return self.idLookup(newkey)
    #fi

    pos = seqPosLookup(dat)
    if pos is None:
      return None
    #fi
    seqid, pos = pos

    cmd = 'curl "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=docsum&log$=seqview" | grep "<title>Homo sapiens" | sed -e \'s/^[[:blank:]]*//g\' | cut -d\  -f4 | cut -f1 -d,' % seqid
    seqnumber = str(utils.getCommandOutput(cmd, shell=True)[0].decode('UTF-8'))
    return (seqnumber.strip(), pos)
  #edef

  def idLookup(self, ID):
    if ID in self.info:
      return self.info[ID]
    else:
      pos = self.__getAssemblyPosition(ID)
      if pos is None:
        return None
      #fi
      self.info[ID] = pos
      return pos
    #fi

  def __call__(self, *pargs, **kwargs):
    return self.lookup(*pargs, **kwargs)
  #edef

  def lookup(self, chromosome, pos, alt, record=False):
    pos = int(pos)
    for r in self.query(chromosome, pos-1, pos):
        if alt in [ ralt.sequence for ralt in r.ALT if hasattr(ralt, "sequence") ]:
            if record:
              return r
            #fi
            if 'RS' in r.INFO:
              return 'rs%d' % r.INFO['RS']
            elif 'SS' in r.INFO:
              return 'ss%d' % r.INFO['SS']
            #fi
        #fi
    #efor
    return None
  #edef

#edef
