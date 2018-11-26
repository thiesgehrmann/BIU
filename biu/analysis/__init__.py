"""
Tools for analysis of different types of data:

 * rnaseq: Analyze RNASeq data
  * Normalization
  * Differential expression (with limma/VOOM)
 * covariates : Analyze covariation between your data and their covariates
"""

from . import rnaseq
from . import covariates
