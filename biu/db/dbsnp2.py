import json

from ..structures import Dataset2
from ..config import settings as settings

from .. import utils
from .. import formats

###############################################################################

class DBSNP2(Dataset2):

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
    
    __assemblySeqIDChromosomes = {}
    
    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        super(DBSNP2, self).__init__("DBSNP/%s" % version , *pargs, **kwargs)
        self.version = version
        
        self._obj.add_file("vcf.vcf.gz", utils.Acquire2().curl(self.versions[self.version]["vcf"]))
        self._obj.add_file("vcf.vcf.gz.tbi", utils.Acquire2().curl(self.versions[self.version]["tbi"]))
        self._obj.register("vcf", ["vcf.vcf.gz", "vcf.vcf.gz.tbi"], lambda f: formats.VCF2(f["vcf.vcf.gz"], tabix=True))
        
        self._obj.add_file("info", utils.Acquire2().touch("info_dict.dbsnp.sqldict.sqlite"))
        self._obj.register("_info", ["info"], lambda f: formats.SQLDict(f["info"]))
    #edef
    
    def filter(self, *pargs, **kwargs):
        return self.vcf.filter(*pargs, **kwargs)
    #edef

    def __getitem__(self, c):
    #Figure out a good way to get info about the rsID from the "pos" file.
    # Tabix doesn't seem to work...
        res = None
        if isinstance(c, int):
            res = self.id_lookup(c)
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
                res = self.id_lookup(c)
            #fi
        else:
            return None
        #fi

        return res
    #edef

    def __getAssemblyPosition(self, ID):
        if not(isinstance(ID, int)):
            utils.msg.error("DBSNP queries must be Integers, not Strings, returning None.")
            return None
        #fi
        utils.dbm("Querying for rs%d via REST." % ID)
        assemblyid = self.versions[self.version]["assemblyid"]

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
        pos = pos + 1 # There seems to be an offset!

        if seqid not in self.__assemblySeqIDChromosomes:
            cmd = 'curl "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=docsum&log$=seqview" | grep "<title>Homo sapiens" | sed -e \'s/^[[:blank:]]*//g\' | cut -d\  -f4 | cut -f1 -d,' % seqid
            seqnumber = str(utils.getCommandOutput(cmd, shell=True).decode('UTF-8'))
            self.__assemblySeqIDChromosomes[seqid] = seqnumber.strip()
        #fi

        return (self.__assemblySeqIDChromosomes[seqid], pos)
    #edef
    
    def id_lookup(self, ID):
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
    #edef

    def __call__(self, *pargs, **kwargs):
        return self.lookup(*pargs, **kwargs)
    #edef

    def lookup(self, chromosome, pos, alt, record=False):
        pos = int(pos)
        for r in self.filter(chromosome, pos-1, pos):
            if alt.lower() in [ a.lower() for a in r.ALT]:
                if record:
                    return r
                #fi
                
                rs = r.INFO.get('RS', None)
                ss = r.INFO.get('SS', None)
                
                if rs is not None:
                    return 'rs%d' % rs
                elif ss is not None:
                    return 'ss%d' % ss
                #fi
            #fi
        #efor
        return None
    #edef

  ###############################################################################

#eclass