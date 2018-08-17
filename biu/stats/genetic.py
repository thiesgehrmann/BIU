from .. import utils
from .. import processing

sstats = utils.py.loadExternalModule("scipy.stats")
np     = utils.py.loadExternalModule("numpy")

from collections import namedtuple

def hardyWeinbergEquilibrium(data, alpha=0.001):
    """
      data : Array of 0/1/2 for allele counts per individual
      alpha: alpha cutoff level
      
      output: HardyWeinbergResult namedtuple
    
    Alpha value taken from:
     https://www.nature.com/articles/ng.2756#method
    """
    
    data = np.array(data)
    
    resTuple = namedtuple("HardyWeinbergResult", [ 'pvalue', 'method', 'chi2statistic', 'significant' ])

    obs    = processing.lst.freq(data)
    pvalue = __SNPHWE(obs.get(1,0), obs.get(0,0), obs.get(2,0))

    #n   = len(data)
    #p   = sum(data) / (sum(data) + sum(data == 0))
    #obs = processing.lst.freq(data)
    #exp = { 0.0 : np.power(1-p,2) * n,
    #        1.0 : 2*p*(1-p) * n,
    #        2.0 : np.power(p,2) * n }

    #chi2 = sum([ np.power(obs.get(c,0) - exp.get(c,0),2)/exp[c] for c in exp ])
    #pvalue = 1 - sstats.chi2.cdf(chi2, 1)
    
    #print(obs, exp, chi2, pvalue)

    #tup = resTuple(pvalue, 'chi2', chi2, pvalue < alpha)
    tup = resTuple(pvalue, 'abecassis_exact', None, pvalue < alpha)
    
    return tup
#edef

def __SNPHWE(obs_hets, obs_hom1, obs_hom2):
    # This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
    # Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
    # Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  

    # NOTE: return code of -1.0 signals an error condition

    # Ported to python by Thies Gehrmann
    
    if (obs_hom1 < 0) | (obs_hom2 < 0) | (obs_hets < 0):
        return -1.0
    #fi

    # total number of genotypes
    N = obs_hom1 + obs_hom2 + obs_hets

    # rare homozygotes, common homozygotes
    obs_homr = min(obs_hom1, obs_hom2)
    obs_homc = max(obs_hom1, obs_hom2)

    # number of rare allele copies
    rare  = obs_homr * 2 + obs_hets

    # Initialize probability array
    # Add a 2 here because I can't be bothered to fix the messy indexing porting from R (where indexing starts at 1)
    probs = np.zeros(2 + rare)

    # Find midpoint of the distribution
    mid = np.floor(rare * ( 2 * N - rare) / (2 * N))
    if ( (mid % 2) != (rare % 2) ):
        mid = mid + 1
    #fi
    mid = int(mid)

    probs[mid + 1] = 1.0
    mysum = 1.0

    # Calculate probablities from midpoint down 
    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr

    while ( curr_hets >=  2):
        probs[curr_hets - 1] = probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
        mysum = mysum + probs[curr_hets - 1]

        # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
        curr_hets = curr_hets - 2
        curr_homr = curr_homr + 1
        curr_homc = curr_homc + 1 
    #ewhile

    # Calculate probabilities from midpoint up
    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr

    while ( curr_hets <= rare - 2):
        probs[curr_hets + 3] = probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum = mysum + probs[curr_hets + 3]

        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        curr_hets = curr_hets + 2
        curr_homr = curr_homr - 1
        curr_homc = curr_homc - 1
    #ewhile

    # P-value calculation
    target = probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    p = min(1.0, sum(probs[probs <= target])/ mysum)
    return p
#edef
