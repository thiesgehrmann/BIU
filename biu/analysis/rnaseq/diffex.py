from ... import utils
from ... import processing
from ... import settings as settings

np = utils.py.loadExternalModule("numpy")
pd = utils.py.loadExternalModule("pandas")

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) 

def voom(formula, expr, covariates, group, contrasts, out_dir=dir_path):
    """
    A wrapper for R differential expression
    Inputs
        formula:    String. An R-style formula
        expr:       Pandas DataFrame of Expression data (columns are genes, rows are samples)
                    Index of dataframe must correspond to index of covariates
        covariates: Pandas Dataframe of Covariates (columns are covariates, rows are samples)
                    Index of dataframe must correspond to index of expression
        group:      The group in which we want to test contrasts (name of column in covariates DataFrame)
        contrasts:  A dictionary of contrasts.
                    e.g. { "c1" : "g2 - g1",
                           "c2" : "g4 - g3",
                           "c3" : "(g4-g3) - (g2-g1)"}

                    Note. The way in which you define the contrasts changes the interpretation!
                    If you say: A : X - Y, and find (from `significant_in`) a gene G is upregulated,
                     it means that G is more highly expressed in X than in Y.

        out_dir:    String. A directory in which files should be output

    Output:
        A pandas dataframe of diff. ex. results

        A typical pipeline would be:

        diffex    = biu.analysis.rnaseq.diffex.voom("~ 0+ condition + covar", expr, covars, "condition", { "A" : "cond2 - cond1"})
        ax        = biu.analysis.rnaseq.diffex.volcanoPlot(diffex)
        sig_genes = biu.analysis.rnaseq.diffex.significant_in(diffex)
        kegg      = biu.db.KEGG()
        enriched  = kegg.enrich(sig_genes['A'])
        
    """
    
    #if not biu.utils.exe.exists('Rscript'):
    #    biu.utils.msg.error('R (and Rscript) is not installed!')
    #    return None
    #fi
    
    Efile = '%s/diffex_voom.expr.csv' % out_dir
    Cfile = '%s/diffex_voom.cov.csv' % out_dir
    Ofile = '%s/diffex_voom.out.csv' % out_dir
    
    print(Efile)
    print(Cfile)
    print(Ofile)
    
    contrast_groups = list(covariates[group].drop_duplicates())
    renamed_contrasts = {}
    for c_name, c_formula in contrasts.items():
        for cg in contrast_groups:
            c_formula = c_formula.replace(cg, 'group_%s' % cg )
        #efor
        renamed_contrasts[c_name] = c_formula
    #efor
    
    contrasts = ';'.join([ "%s=%s" % (c_name, c) for (c_name, c) in renamed_contrasts.items() ])
    print(contrasts)
    
    # Change all numerical covariates to strings (Otherwise R reads them as numeric covariates)
    C = covariates.copy()
    for column in C:
        if column == group:
            C[column] = C[column].apply(lambda x: 'group_%s' % str(x))
        elif C[column].dtype.name  == 'category':
            C[column] = C[column].apply(lambda x: 'cat_%s' % str(x))
        #fi
    #efor
    samples_renamed = { s : 'sample_%d' % (i+1) for (i,s) in enumerate(C.index) }
    covars_renamed  = { c : 'covar_%d' % (i+1) for (i,c) in enumerate(C.columns) }
    C = C.rename(columns=covars_renamed, index=samples_renamed) 
    group = covars_renamed[group]

    for (c, cr) in covars_renamed.items():
        formula = formula.replace(c, cr)
    #efor
    
    print(formula)
    
    # Implement some kind of check that the parameters are in the covariates file...
    
    C.to_csv(Cfile, sep=';', index_label='index')

    E = expr.copy()
    gene_renamed = { g: 'gene_%d' % (i+1) for (i,g) in enumerate(E.columns) }
    gene_denamed = { v:k for (k,v) in gene_renamed.items() }
     
    E = E.rename(columns=gene_renamed, index=samples_renamed)
    E.to_csv(Efile, sep=';', index_label='index')    
    
    cmd = 'Rscript "%s" "%s" "%s" "%s" "%s" "%s" "%s"' % ('%s/diffex_voom.R' % dir_path, formula, group, contrasts, Efile, Cfile, Ofile)
    print(cmd)
    p = utils.exe.runCommand(cmd)
    if p != 0:
        utils.msg.error('R returned an error. Please check input.')
        return None
    #fi
    
    R = pd.read_csv(Ofile, sep=';').rename(columns={'adj.P.Val' : 'qvalue', 'P.Value' : 'pvalue'})
    R['gene'] = R.gene.apply(lambda g: gene_denamed[g])
    return R
#edef

def summary(testRes, alpha=0.05, lfcThresh=0.5, col_contr='contr', col_pval='qvalue', col_lfc='logFC'):
    """
    summary: Returns a summary of the differential expression tests
    Inputs: 
        testRes : The table of tests (e.g. from voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
            
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output: List of genes significant by the specified conditions
    """
    R = testRes.copy()
    R['significant'] = 1 * (R[col_pval].values < alpha) & (np.abs(R[col_lfc].values) > lfcThresh )
    R['up']          = 1 * ((R.significant == 1) & (R[col_lfc] > 0))
    R['down']        = 1 * ((R.significant == 1) & (R[col_lfc] < 0))
    
    return R.groupby(col_contr).agg(sum)[['significant','up','down']]
#edef


def sigTests(testRes, alpha=0.05, lfcThresh=0.5, col_pval='qvalue', col_lfc='logFC', col_contr='contr'):
    """
    Return the set rows of significant tests.
    Input:
        testRes: the output of a differential expression test (e.g. voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        A subset of the significant rows in input testRes
    """
    return testRes[(testRes[col_pval] < alpha) & (testRes[col_lfc].abs() > lfcThresh)]
#edef

def significant_in(diffex, contr=None, split_updown=False, alpha=0.05, lfcThresh=0.5,
                   col_contr='contr', col_pval='qvalue', col_lfc='logFC', col_index='gene'):
    """
    Return a dictionary of significant genes per contrast
    Inputs:
        diffex: The output of a differential expression (e.g. voom)
        contr: String or list of strings of contrast names (or None)
        split_updown: Split the differential expressions by up/down expressed
    Outputs:
        dict of { contr : genes } if split_updown is False
        dict of { contr: { 'up' : genes, 'down' : genes } } if split_updown is True
    """
    if contr is None:
        contr = list(diffex[col_contr].drop_duplicates())
    #fi
    if isinstance(contr, str):
        contr = [ contr ]
    #fi
    
    sigs = sigTests(diffex,alpha=alpha, lfcThresh=lfcThresh,
                    col_contr=col_contr, col_lfc=col_lfc, col_pval=col_pval)
    sigin = {}
    if split_updown:
        sigin = { (c,direction): sigs[(sigs[col_contr] == c) & ((sigs[col_lfc] > 0) if direction == 'up' else (sigs[col_lfc] <= 0) )][col_index].values
                  for c in contr for direction in ['up', 'down']}
    else:
        sigin = { c : sigs[sigs[col_contr] == c][col_index].values for c in contr }
    #fi
    return sigin
#edef

def isSigIn(testRes, contr=None, logic=all, alpha=0.05, lfcThresh=0.5,
                  col_contr='contr', col_index='gene', col_pval='qvalue', col_lfc='logFC'):
    """
    isSigIn: Returns genes that are significant in a given set of contrasts.
    Inputs: 
        testRes : The table of tests
        contr: [list(str)] The contrans in which we wish to inspect significance
        logic: [Bool -> Bool] How to combine, e.g. any means sigificant in at least one contrast
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
            
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output: List of genes significant by the specified conditions
    """
    
    if contr is None:
        contr = set(testRes[col_contr].values)
    #fi
    relRes = sigTests(testRes, alpha=alpha, lfcThresh=lfcThresh, col_pval=col_pval,
                      col_lfc=col_lfc, col_contr=col_contr)
    sigGeneGroups = processing.lst.group(relRes[[col_contr, col_index]].values, key=lambda x: x[1], value=lambda x: x[0])
    sigGenes = [ g for g in sigGeneGroups if logic([ (c in sigGeneGroups[g]) for c in contr ]) ]
    return sigGenes
#edef

def volcanoPlot(testRes, alpha=0.05, lfcThresh=0.5, contr=None, ax=None,
                col_pval='qvalue', col_lfc='logFC', col_contr='contr', **kwargs):
    """
    volcanoPlot: Plots volcanoplots for each contrast
    Inputs:
        testRes : The table of tests (e.g. output from voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        contr: String or list of Strings. Which contrast to investigate? 
        ax: Matplotlib axes to plot on. Same length as contr
        
        
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        List of genes significant by the specified conditions
    """

    tests = testRes[col_contr].drop_duplicates().values

    if contr is not None:
        tests = [ contr ] if isinstance(contr, str) else contr
    #fi
    
    axes = None
    if ax is not None:
        axes = ax if hasattr(ax, '__len__') else [ ax ]
        assert len(ax) == len(tests)
    else:
        ncols = min(np.ceil(np.sqrt(len(tests))), 4)
        nrows = np.ceil(len(tests) / ncols)
        fig, axes = utils.figure.subplots(ncols=int(ncols), nrows=int(nrows),
                                          figsize = (5*ncols, 5*nrows),
                                          sharey=True, sharex=True, **kwargs)
    #fi

    rangeFC = (testRes[col_lfc].min(), testRes[col_lfc].max())
    rangePV = -np.log10([testRes[col_pval].max(), testRes[col_pval].min()])
    sigThresh = -np.log10(alpha)

    for idx, test in enumerate(tests):
        relRows = testRes[testRes[col_contr] == test]
        sigRows = sigTests(relRows, alpha=alpha, lfcThresh=lfcThresh,
                           col_pval=col_pval, col_lfc=col_lfc, col_contr=col_contr)
        axes[idx].scatter(relRows[col_lfc], -np.log10(relRows[col_pval]), s=0.5)
        axes[idx].set_title(test)
        axes[idx].plot([rangeFC[0], rangeFC[1]], [sigThresh, sigThresh], c='orange')
        axes[idx].plot([-lfcThresh, -lfcThresh], rangePV, c='orange')
        axes[idx].plot([lfcThresh, lfcThresh], rangePV, c='orange')
        axes[idx].scatter(sigRows.logFC, -np.log10(sigRows[col_pval]), c='r', s=0.5)
        axes[idx].set_ylabel('-log pvalue')
        axes[idx].set_xlabel('log Fold Change')
    #efor
    
    return fig, axes
#edef

def compair(diffex, contrA, contrB, alpha=0.05, lfcThresh=0.5,
            col_contr='contr', col_index='gene', col_pval='qvalue', col_lfc='logFC'):
    """
    Compare a pair (compair) of contrasts, per gene
    Inputs:
        diffex: The table of tests (e.g. output from voom)
        contrA: The baseline contrast
        contrB: The contrast to compare to
        alpha: pvalue threshold
        lfcThresh: LogFC threshold        
        
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
        
    Outputs:
        Dataframe of comparisons, per gene
        
        Comparisons are given:
            significant: Is the gene significantly differentially expressed in any condition (pvalue < alpha) & lfc > lfcThresh
            consistent:  Does the Log Fold change have the same sign in both contrasts
            stronger:    Is the Log Fold Change consistent and larger in contrast B than in contrast A?
            weaker:      Is the LFC consistent and smaller in contrast B than in contrast A?
            slope:       In terms of the volcano plot, with which slope does the gene move?
            direction:   In terms of LFC, which direction does the gene move? (-1 to left, +1 to right)
    """
    
    D = diffex[diffex[col_contr].isin([contrA, contrB])].copy()
    D = D.pivot(index='gene', columns='contr', values=['qvalue', 'logFC'])
    D['significant'] = list(map(np.any, (D[col_pval].values < 0.05) & (D[col_lfc].abs().values > lfcThresh )))
    D['consistent']  = list(map(lambda x: ~np.logical_xor(*x), (D[col_lfc] > 0).values))
    D['stronger']    = D.consistent.values & (D[col_lfc][contrB].abs().values > D[col_lfc][contrA].abs().values)
    D['weaker']      = D.consistent.values & (D[col_lfc][contrB].abs().values < D[col_lfc][contrA].abs().values)
    D['slope']       = (-np.log10(D[col_pval][contrB].values)+np.log10(D[col_pval][contrA].values))/(D[col_lfc][contrB].values-D[col_lfc][contrA].values)
    D['direction']   = np.sign(D[col_lfc][contrB].values - D[col_lfc][contrA].values)
    
    return D
#edef

def pairedVolcanoPlot(diffex, contrA, contrB, alpha=0.05, lfcThresh=0.5, only_significant=True, ax=None,
                      color='stronger', col_contr='contr', col_index='gene', col_pval='qvalue', col_lfc='logFC'):
    """
    pairedVolcanoPlot: Plot two volcanoplots on top of each other, and join genes by lines
    Inputs:
        diffex: The table of tests (e.g. output from voom)
        contrA: Contrast A to display (in blue)
        contrB: Contrast B to display (in red)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        only_significant: Only plot genes that are significant in either one of the two contrasts
        color: How to color the plots between the two contrasts?
               See output from compair().
        ax: Matplotlib axes to plot on. Same length as contr
        
        
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        List of genes significant by the specified conditions
    """
    
    cmp = compair(diffex, contrA, contrB, alpha=alpha, lfcThresh=lfcThresh,
                  col_contr=col_contr, col_index=col_index, col_pval=col_pval, col_lfc=col_lfc)
    
    if only_significant:
        cmp = cmp[cmp.significant]
    #fi
        
    if ax is None:
        fig, axes = utils.figure.subplots(ncols=1, nrows=1)
        ax = axes[0]
    #fi
    
    ax.scatter(cmp[col_lfc][contrA].values, -np.log10(cmp[col_pval][contrA].values), c='blue', label=contrA, s=0.5)
    ax.scatter(cmp[col_lfc][contrB].values, -np.log10(cmp[col_pval][contrB].values), c='red',  label=contrB, s=0.5)
    
    lines = ax.plot([cmp[col_lfc][contrA].values, cmp[col_lfc][contrB].values],
                    -np.log10([cmp[col_pval][contrA].values, cmp[col_pval][contrB].values]),
                    alpha=0.5, linewidth=0.1)
    

    for l, c in zip(lines, cmp[color]):
        l.set_color('g' if 1*c > 0 else 'r')
    #efor
    
    rangeFC = (np.min(cmp[col_lfc]),np.max(cmp[col_lfc]))
    rangePV = -np.log10((np.min(cmp[col_pval]),np.max(cmp[col_pval])))
    sigThresh = -np.log10(alpha)
        
    ax.plot([rangeFC[0], rangeFC[1]], [sigThresh, sigThresh], c='orange', alpha=0.5, linewidth=1)
    ax.plot([-lfcThresh, -lfcThresh], rangePV, c='orange', alpha=0.5, linewidth=1)
    ax.plot([lfcThresh, lfcThresh], rangePV, c='orange', alpha=0.5, linewidth=1)
    ax.legend()
    
    ax.set_title('Paired volcano plot\n%s (blue) vs %s (red)' % (contrA, contrB))
    ax.set_xlabel('Log Fold Change')
    ax.set_ylabel('- log10 pvalue')
    
    return ax, cmp
#edef

