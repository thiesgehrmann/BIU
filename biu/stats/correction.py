def bonferroni(pvals):
  """
    Bonferroni FWER correction of p-values
    Inputs:
     - pvals: An array of p-values
    Outputs:
     - Corrected q-values
  """
  return np.array(pvals) * len(pvals)
#edef

def fdr(p, alphai=0.05, type='bh'):

  nt     = len(p);  # Number of tests
  ps   = sorted(p);
  indx = [ i[0] for i in sorted(enumerate(p), key=lambda x:x[1]) ];

  if type == 'bhy':
    cm    = sum([ 1.0/float(i) for i in xrange(nt)] );
    klist = [ (float(i+1)/(float(nt) * cm)) for i in xrange(nt) ];
  else:
    klist = [ (float(i+1)/float(nt)) for i in xrange(nt) ];
  #fi

    # Adjust pvalues to qvalues
  padj = [ ps[i] / klist[i] for i in xrange(nt)];
    # Fix pvalues larger than 1
  q = [ qi if qi < 1.0 else 1.0 for qi in padj ];

    # Monotonicity
  qm = [];
  prev_v = q[0];
  for v in q:
    qm.append(max(prev_v, v));
    prev_v = qm[-1];
  #efor

    # get back to original sorting
  qrs = [0] * nt;
  for i in xrange(nt):
    qrs[indx[i]] = qm[i];
  #efor

  return qrs;
#edef


