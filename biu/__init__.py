"""
BIU: BioInformatics Utilities.

BIU is a library of utilities and tools that provides:

  * Access to several standard (and not so standard) Bioinformatics formats (`biu.formats`)

  * Automated retrieval and interface to data from several databases (`biu.db`):
    * Genomes (Ensembl, Flybase, Wormbase)
    * Database/Datasets (ClinVar, DBSnp, GnomAD, Kegg, UniProt, CADD, etc.)

  * Mapping between genome Identifiers (`biu.maps`):
    * Within a genome, where possible (gene -> Transcript -> Protein)
    * Between databases (with Ensembl BioMart)
    * Between species (DIOPT)

  * Bioinformatic pipelines (`biu.pipelines`):
    * BLAST/Diamond (`biu.pipelines.Blast`)
    * Variant Effect Predictor (biu.pipelines.VEP)
    * Write your own! (`biu.structures.Pipeline`)

  * Analysis of bioinformatic data (`biu.analysis`)
    * RNASeq data: (`biu.analysis.rnaseq`)
    * Covariates: (`biu.analysis.covariates`)

  * Medical tools (`biu.medical`):
    * health parameters: (`biu.medical.health`) (Framingham Risk Score, etc)

  * Data Structure processing (`biu.processing`):
    * Python lists
    * numpy Matrices
    * pandas DataFrames
   
  * Statistical analysis (`biu.stats`)
    * Regression (`biu.stats.regression`):
      * LOWESS regression
    * Enrichment (`biu.stats.enrichment`)
    * Multiple Testing Correction (FWER, FDR) (`biu.stats.correction`)
    * Statistical genetics (`biu.stats.genetics`)
      * Hardy Weinberg Equilibrium test)

  * Mathematical functionality (`biu.math`)
    * Interpolation (`biu.math.interpolation`)
    

BIU depends on a lot of modules.
However, much functionality in BIU does not depend on many other packages.
Additionally, BIU is implemented in such a way that it will load, even if none of those dependencies are actually installed.
This means that you can load BIU without any dependencies installed, and check if the functionality you are interested in works without them.
A recommended baseline for most functionality would be:
  * numpy
  * pandas
  * scipy
  * tabix (needed for some datasets, and VCF indexing)
  * vcf (needed for some datasets and VCF loading)
  * intervaltree (needed for indexing formats such as VCF and GFF files)

All external dependencies currently used are:
  * fastcluster
  * intervaltree
  * matplotlib
  * matplotlib_venn
  * numpy
  * openpyxl
  * pandas
  * scipy
  * seaborn
  * sklearn
  * snakemake
  * sqlite3
  * statsmodels
  * tabix
  * vcf
  * xlrd

"""

# Make settings available
from .config import config
from .config import settings

# Make internal structures and utilities available
from . import structures as structures
from . import utils as utils

# Make formats accesible
from . import formats as formats

# Make genomes, datasets and Identifier mappings acessible
from . import genomes as genomes
from . import db as db
from . import maps as maps

# Make higher-level mathematical and statistical operations available
from . import math as math
from . import stats as stats

# Make data structure processing available
from . import processing as processing

# Make tools for analysis and pipelines available
from . import analysis as analysis
from . import pipelines as pipelines
from . import medical as medical

def __version__():
  print("BIU (Bio Utilities) python module")
  print(config.settings.dumps())
  print(" Current config hash: %s" % processing.lst.hash(config.settings.dumps()))
#edef

__missingDependencies = settings.missingDependencies()
if len(__missingDependencies[0]) > 0:
  utils.msg.warning("The following dependencies of BIU are missing. Functionality of BIU will be affected.\n  %s" % ', '.join(__missingDependencies[0]))
#fi

if len(__missingDependencies[1]) > 0:
  utils.msg.warning("Some optional dependencies of BIU are missing. Functionality of BIU may be affected.\n  %s" % ', '.join(__missingDependencies[1]))
#fi
