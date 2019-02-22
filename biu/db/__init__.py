## Old Dataset databases
from .caddUtils import CADD as CADD

from .clinVarUtils import ClinVar as ClinVar

from .gnomadUtils import Gnomad as Gnomad

from .gtexUtils import GTeX as GTeX

from .uniprotUtils import UniProt as UniProt

from .hagrUtils import HAGR as HAGR

from .llsUtils import LLS as LLS
from .bbmriUtils import BBMRI as BBMRI

from .goto import GOTO as GOTO

from .goUtils import GO as GO

from .keggUtils import KEGG as KEGG

from .dbsnpUtils import DBSNP as DBSNP

### Dataset2 databases

from .iris import Iris
from .gwas_catalog import GWAS_Catalog
from .cosmic import Cosmic as Cosmic
from .lls2 import LLS2 as LLS2 # A re-implementation of LLS
from .bbmri2 import BBMRI2 as BBMRI2 # A re-implementation of BBMRI
from .gnomad2 import Gnomad2 as Gnomad2


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

__datasets = [ GO, LLS, GOTO, Cosmic, HAGR, KEGG, DBSNP, BBMRI, CADD, ClinVar, UniProt, Gnomad, GWAS_Catalog, Iris, GTeX ]

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


