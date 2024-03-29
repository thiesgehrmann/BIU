"""
Tools for analysis of different types of data:

 * rnaseq: Analyze RNASeq data
  * Normalization
  * Differential expression (with limma/VOOM)
 * covariates : Analyze covariation between your data and their covariates
 * trajectory: Analyze the movement of samples through PC space
"""

from . import hierarchy
from . import rnaseq
from . import covariates
from . import trajectory
from . import metabolomics
from .trajectory import Trajectory
from . import microbiome
from . import clustering
