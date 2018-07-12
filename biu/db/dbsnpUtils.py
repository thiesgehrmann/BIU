
import json

from ..structures import Dataset
from ..config import settings as settings

from .. import utils
from .. import formats

###############################################################################

class DBSNP(Dataset):

  version = None
  where    = None
  fileIndex = None

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

  def __init__(self, version=list(versions.keys())[0], where=None):
    fileIndex = self.__genFileIndex(version)
    Dataset.__init__(self, fileIndex)
    self.version = version

    self._addStrFunction( lambda s: "Version : %s" % s.version )
    self._registerObject("_vcf", formats.VCF, [ "vcf", "tbi" ], fileIndex["vcf"].path, tabix=True, )
    self._registerObject("_info", formats.SQLDict, [ "info"], fileIndex["info"].path)

  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/dbSNP/%s' % ( (settings.getDataDir() if where is None else where), version)
     vData = self.versions[version]
     files = {}
     files['vcf'] = utils.Acquire(where=where).curl(vData["vcf"]).finalize('%s/all_snps.vcf.bgz' % finalPath)
     files['tbi'] = utils.Acquire(where=where).curl(vData["tbi"]).finalize('%s/all_snps.vcf.bgz.tbi' % finalPath)
     files['info'] = utils.Acquire(where=where).touch("%s/info_dict.sqldict.sqlite" % finalPath)
     return files
  #edef

  ###############################################################################

  def query(self, *pargs, **kwargs):
    return self._vcf.query(*pargs, **kwargs)
  #edef

  def queryRegions(self, *pargs, **kwargs):
    return self._vcf.queryRegions(*pargs, **kwargs)
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

  __assemblySeqIDChromosomes = {}

  def __getAssemblyPosition(self, ID):
    if not(isinstance(ID, int)):
      utils.msg.error("DBSNP queries must be Integers, not Strings, returning None.")
      return None
    #fi
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
      utils.warning("Cannot find rs%d for genome %s!" % (ID, assemblyid))
      return None
    #edef

    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/%d' % ID
    dat = str(utils.getCommandOutput('curl -X GET --header "Accept: application/json" "%s"' % url).decode('UTF-8'))
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

    if seqid not in self.__assemblySeqIDChromosomes:
      cmd = 'curl "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=docsum&log$=seqview" | grep "<title>Homo sapiens" | sed -e \'s/^[[:blank:]]*//g\' | cut -d\  -f4 | cut -f1 -d,' % seqid
      seqnumber = str(utils.getCommandOutput(cmd, shell=True).decode('UTF-8'))
      self.__assemblySeqIDChromosomes[seqid] = seqnumber.strip()
    #fi

    return (self.__assemblySeqIDChromosomes[seqid], pos)
  #edef

  def idLookup(self, ID):
    if ID in self._info:
      res = self._info[ID]
      if res == -1:
        return (None, None)
      #fi
      return res
      
    else:
      pos = self.__getAssemblyPosition(ID)
      if pos is None:
        utils.dbm("Pos is None!")
        self._info[ID] = -1
        return (None, None)
      #fi
      self._info[ID] = pos
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
