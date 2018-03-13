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
from .cosmicUtils import listVersions as listCosmicVersions

from .hagrUtils import HAGR as HAGR
from .hagrUtils import listVersions as listHagrVersions

from .reactomeUtils import Reactome as Reactome
from .reactomeUtils import listVersions as listReactomeVersions

from .neo4j import Neo4jDB as Neo4jDB
from .neo4j import listVersions as listNeo4jVersions

from .llsUtils import LLS as LLS
from .llsUtils import listVersions as listLLSVersions

from .bbmriUtils import BBMRI as BBMRI
from .bbmriUtils import listVersions as listBBMRIVersions

from .goUtils import GO as GO
from .goUtils import listVersions as listGOVersions

from .keggUtils import KEGG as KEGG
from .keggUtils import listVersions as listKEGGVersions

from .mirmineUtils import MiRmine as MiRmine
from .mirmineUtils import listVersions as listMiRmineVersions

def list():
  print("Available databases:")
  for db in sorted(["Genomes", "CADD", "ClinVar", "Gnomad", "GTeX", "UniProt", "Cosmic", "HAGR", "Reactome", "LLS", "BBMRI", "GO", "KEGG", "MiRmine" ]):
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

  print("\nReactome:")
  listReactomeVersions()

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
