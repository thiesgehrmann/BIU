from .. import utils

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

###############################################################################

def correct(pvals, correctionType, **kwargs):
  """
  correct: Multiple testing correction.

  Inputs:
   - pvals: list of p values
   - correctionType: Type of correction
     - fwer|bonferroni => bonferroni()
     - fdr|fdr_bh => fdrBH()
     - fdr_bhy => fdrBHY()
   - **kwargs: Additional arguments (See arguments to called functions
 
  Outputs:
   - corrected p-values
  """
  correctionType = correctionType.lower()

  if correctionType in [ 'fwer', 'bonferroni' ]:
    q = bonferroni(pvals)
  elif correctionType in [ 'fdr', 'fdr_bh' ]:
    q = fdrBH(pvals)
  elif correctionType in [ 'fdr_bhy' ]:
    q = fdrBHY(pvals, **kwargs)
  else:
    utils.warning("correctionType='%s' is unknown. Falling back to bonferroni." % str(correctionType))
    q = bonferroni(pvals)
  #fi

  return q
#edef

def bonferroni(pvals):
  """
    Bonferroni FWER correction of p-values
    Inputs:
     - pvals: An array of p-values
    Outputs:
     - Corrected q-values
  """
  q = np.array(pvals) * len(pvals)
  q[q > 1] = 1
  return q
#edef

def fdrBH(pvals):
  """
    Benjamini-Hochberg correction of p-values
    Inputs:
     - pvals: An array of p-values
    Outputs:
     - Corrected q-values
  """
  nt = len(pvals)
  q  = np.zeros(nt)
  sp = sorted(enumerate(pvals), key=lambda x: x[1])
  prevq = 0.0
  for i, (idx, pv) in enumerate(sp):
    # Correct p-value, and ensure monotonicity
    pcorr  = pv * (nt/(i+1))
    q[idx] = max(pcorr, prevq)
    prevq = q[idx]
  #efor
  q[q > 1] = 1
  return q
#edef

def fdrBHY(pvals, cm=None):
  """
    Benjamini-Hochberg-Yekutieli correction of p-values
    Inputs:
     - pvals: An array of p-values
     - cm : Independence factor np.sum(1.0/(np.array( range(nt)) + 1 )), nt=len(pvals)
    Outputs:
     - Corrected q-values
  """
  nt = len(pvals)
  q  = np.zeros(nt)
  sp = sorted(enumerate(pvals), key=lambda x: x[1])
  prevq = 0.0
  if cm is None:
    cm = np.sum(1.0/(np.array( range(nt)) + 1 ));
  #fi
  for i, (idx, pv) in enumerate(sp):
    # Correct p-value, and ensure monotonicity
    pcorr  = pv * ((nt*cm)/(i+1))
    q[idx] = max(pcorr, prevq)
    prevq = q[idx]
  #efor
  q[q > 1] = 1
  return q
#edef

