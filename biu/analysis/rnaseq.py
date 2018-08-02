from .. import utils

sstats = utils.py.loadExternalModule("scipy.stats")
np     = utils.py.loadExternalModule("numpy")
pd     = utils.py.loadExternalModule("pandas")

def normalize(method, E, *pargs, **kwargs):
    methods = { 'voom' : __normalize_voom,
                'tmm' : __normalize_tmm,
                'fpkm' : __normalize_fpkm,
                'nonzero' : lambda x: x}
    method = method.lower()
    if method in methods:
        D = E.values
        notNullCols = np.where(np.any(D > 0, axis=0))[0]
        notNullRows = np.where(np.any(D > 0, axis=1))[0]
        utils.msg.dbm("Nonzero rows, columns: %d, %d" % (len(notNullRows), len(notNullCols)))
        D = D[notNullRows,:][:,notNullCols]
        N = methods[method](D, *pargs, **kwargs)
        
        return pd.DataFrame(N, columns=E.columns[notNullCols], index=E.index[notNullRows])
    else:
        raise NotImplementedError('%s normalization is not implemented' % method)
    #fi
#edef

def __normalize_voom(E):
    """log2 CPM normalization
       log2CPM = log_2 (10^6 * (r + 0.5) / (R + 1.0)) = 6 log_2 10 + log_2 (r+.5) - log_2 (R + 1)"""
    t = np.log2(E + 0.5) + 6 * np.log2(10) - np.log2(E.sum(axis=1) + 1.0)[:,np.newaxis]
    return t
    
def __normalize_tmm(E, refID=None, logratioTrim=0.3, sumTrim=0.05, aCutoff=-1e10):
    """ TMM normalization
        As implemented: https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
        However, here implemented in python
    """
    if refID is None:
        uq  = np.apply_along_axis(lambda x: np.percentile(x, 75), axis=1, arr=E)
        muq = np.mean(uq)
        refID, _ = min(enumerate(np.abs(uq-muq)), key=lambda x: x[1])
    #fi
    
    def calcNormFactors(sampleID):
        sample = E[sampleID, :]
        ref    = E[refID, :]
        
        nS = sum(sample)
        nR = sum(ref)

        logR = np.log2((sample/nS) / (ref/nR))
        absE = (np.log2(sample/nS) + np.log2(ref/nR))/2
        v    = (nS-sample)/(nS*sample) + (nR-sample)/(nR*sample)
        
        fin = np.isfinite(logR) & np.isfinite(v) & (absE > aCutoff)
        
        logR = logR[fin]
        absE = absE[fin]
        v    = v[fin]
        
        if np.max(np.abs(logR)) < 1e-6:
            return 0
        #fi
        
        #taken from the original mean() function
        n = len(logR)
        loL = np.floor(n * logratioTrim) + 1
        hiL = n + 1 - loL
        loS = np.floor(n * sumTrim) + 1
        hiS = n + 1 - loS
        #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
        #a fix from leonardo ivan almonacid cardenas, since rank() can return
        #non-integer values when there are a lot of ties
        
        logRank = np.array(sstats.rankdata(logR))
        absRank = np.array(sstats.rankdata(absE))
        keep = (logRank>=loL) & (logRank<=hiL) & (absRank>=loS) & (absRank<=hiS)

        fNumin = logR[keep]/v[keep]
        fDenom = 1/v[keep]
        f = sum(fNumin[np.isfinite(fNumin)]) / sum(fDenom[np.isfinite(fDenom)])
        
        #Results will be missing if the two libraries share no features with positive counts
        #In this case, return unity
        if not np.isfinite(f):
            f = 0
        #efor
        return f
    #edef
    
    f = np.array([ calcNormFactors(k) for k in range(E.shape[0])])
    f = np.power(2, f)
    
    #f = f/np.exp(np.mean(np.log(f)))
    utils.msg.dbm("TMM NormFactors (Ref=%d) = %s" % (refID, str(f)))
    return E / f[:,np.newaxis]
#edef
    
def __normalize_fpkm(E, gff):
    raise NotImplementedError("FPKM normalization isn't implemented yet.")
#edef
    

