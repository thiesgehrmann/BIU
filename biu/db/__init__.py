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

def list():
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

#edef
