from .caddUtils import CADD as CADD
from .caddUtils import listVersions as listCADDVersions

from .clinVarUtils import ClinVar as ClinVar
from .clinVarUtils import listVersions as listClinVarVersions

from .gnomadUtils import Gnomad as Gnomad
from .gnomadUtils import listVersions as listGnomadVersions

from .genomeUtils import Genome as Genome
from .genomeUtils import listGenomes as listGenomes

from .gtexUtils import GTeX as GTeX
from .gtexUtils import listVersions as listGTeXVersions

from .uniprotUtils import UniProt as UniProt
from .uniprotUtils import listVersions as listUniProtVersions

from .cosmicUtils import Cosmic as Cosmic

from .hagrUtils import HAGR as HAGR

from .reactomeUtils import Reactome as Reactome
from .reactomeUtils import listVersions as listReactomeVersions

#from .neo4j import Neo4jDB as Neo4jDB
#from .neo4j import listVersions as listNeo4jVersions

from .llsUtils import LLS as LLS

from .bbmriUtils import BBMRI as BBMRI
from .bbmriUtils import listVersions as listBBMRIVersions

from .goUtils import GO as GO

from .keggUtils import KEGG as KEGG
from .keggUtils import listVersions as listKEGGVersions

from .mirmineUtils import MiRmine as MiRmine
from .mirmineUtils import listVersions as listMiRmineVersions

from .dbsnpUtils import DBSNP as DBSNP
from .dbsnpUtils import listVersions as listDBPSNPVersions

from .rvisUtils import RVIS as RVIS
from .gdiUtils import GDI as GDI
from .dneUtils import DNE as DNE

__datasets = [ GO, LLS, Cosmic, HAGR ]

def versions(db = None):
  if db is None:
    db = __datasets
  else:
    db = [ db ]
  #fi

  for d in sorted(db, key=lambda x: x.__name__):
    print(d.__name__)
    for version in d.versions:
      print(" * %s" % version)
    #efor
  #efor
#edef

def list():
  print("Available databases:")
  for db in sorted(["Genomes", "CADD", "ClinVar", "Gnomad", "GTeX", "UniProt", "Cosmic", "HAGR", "LLS", "BBMRI", "GO", "KEGG", "MiRmine", "RVIS", "GDI", "DNE" ]):
    print(" * %s" % db)
  #efor
#eclass

def listVersions():
  print("CADD:")
  listCADDVersions()

  print("\nClinVar:")
  listClinVarVersions()

  print("\nGnomAD")
  listGnomadVersions()

  print("\nGenomes:")
  listGenomes()

  print("\nGTeX:")
  listGTeXVersions()

  print("\nUniProt:")
  listUniProtVersions()

#  print("\nReactome:")
#  listReactomeVersions()

  print("\nNeo4j:")
  listNeo4jVersions()

  print("LLS:")
  listLLSVersions()

  print("BBMIR:")
  listBBMRIVersions()

  print("GO:")
  listGOVersions()

  print("KEGG:")
  listKEGGVersions()

  print("MiRmine:")
  listMiRmineVersions()

#edef
