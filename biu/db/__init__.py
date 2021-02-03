## Old Dataset databases, existing only for backwards compatability

from .bbmriUtils import BBMRI as BBMRI
from .caddUtils import CADD as CADD
from .gnomadUtils import Gnomad as Gnomad
from .goUtils import GO as GO
from .goto import GOTO as GOTO
from .llsUtils import LLS as LLS
from .dbsnpUtils import DBSNP as DBSNP


## Old dataset databases, not yet existing in dataset2 format
from .clinVarUtils import ClinVar as ClinVar
from .gtexUtils import GTeX as GTeX
from .hagrUtils import HAGR as HAGR
from .keggUtils import KEGG as KEGG
from .uniprotUtils import UniProt as UniProt

### Dataset2 databases

from .bbmri2 import BBMRI2 as BBMRI2 # A re-implementation of BBMRI
from .cadd2 import CADD2 as CADD2
from .cosmic import Cosmic as Cosmic
from .genage import GenAge as GenAge
from .dbsnp2 import DBSNP2 as DBSNP2 # a re-implementation of DBSNP
from .disgenet import DisGeNet as DisGeNet
from .gnomad2 import Gnomad2 as Gnomad2
from .go2 import GO2 as GO2
from .goto2 import GOTO2 as GOTO2
from .gwas_catalog import GWAS_Catalog
from .hgnc import HGNC
from .iris import Iris
from .lls2 import LLS2 as LLS2 # A re-implementation of LLS
from .ncbi_taxonomy import NCBITaxonomy as NCBITaxonomy
from .reactome import Reactome as Reactome
from .stringdb import StringDB as StringDB
from .tfcat import TFCAT as TFCAT
from .wormbase import WormBase as WormBase




## SUPER OLD ResourceManager/LazyObject databases

#from .rvisUtils import RVIS as RVIS
#from .gdiUtils import GDI as GDI
#from .dneUtils import DNE as DNE
#from .mirmineUtils import MiRmine as MiRmine
#from .mirmineUtils import listVersions as listMiRmineVersions
#from .reactomeUtils import Reactome as Reactome
#from .reactomeUtils import listVersions as listReactomeVersions

#from .neo4j import Neo4jDB as Neo4jDB
#from .neo4j import listVersions as listNeo4jVersions
#from .genomeUtils import Genome as Genome
#from .genomeUtils import listGenomes as listGenomes


__datasets = [ BBMRI, CADD, Gnomad, GO, GOTO, LLS, # Old format, also in new format
               ClinVar, DBSNP, GTeX, HAGR, KEGG, UniProt, # Old format, not yet in new format
               BBMRI2, CADD2, Cosmic, Gnomad2, GO2, GOTO2, GWAS_Catalog, HGNC,
               Iris, LLS2, NCBITaxonomy, Reactome, StringDB, TFCAT  ] # New format
               

def versions(db = None):
    if db is None:
        db = __datasets
    else:
        db = [ db ]
    #fi

    for d in sorted(db, key=lambda x: x.__name__):
        print(d.__name__)
        if hasattr(d, 'versions'):
            for version in d.versions:
                print(" * %s" % version)
            #efor
        #fi
    #efor
#edef

def list():
    print("Available databases:")
    for db in sorted(__datasets):
        print(" * %s" % db.__name__)
#efor


