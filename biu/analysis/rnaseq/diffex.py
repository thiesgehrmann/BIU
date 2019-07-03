from ... import utils
from ... import ops
from ... import settings as settings
from ... import stats
from ... import R

np = utils.py.loadExternalModule("numpy")
pd = utils.py.loadExternalModule("pandas")

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

_defaults = dict(col_contr='contr',
                 col_index='gene',
                 col_pval='pvalue',
                 col_qval='fdr',
                 col_lfc='logFC',
                 lfcThresh=0.5,
                 alpha=0.05)

###################################################################################

def limma(formula, cov, expr, contrasts=None, random_effect=None, voom=True):
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
        random_effect: String. The name of a column to define a random effect (e.g. per person).
        voom: Boolean. Use voom normalization (Default true)
        
    returns:
    A dataframe
    """
    r = R()
    r.push(cov=cov, expr=expr, f=formula, re=random_effect, voom=voom)
    r.exec('source(biu$limma.tools')
    r.exec('diffex <- limma.rnaseq(formula(f), cov, expr, contrasts, re, voom)')
    return r.get('diffex')
#edef

###################################################################################

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
        NOTE: The qvalue column is the corrected p-value provided by limma. This is corrected WITHIN each condition
              The fdr column gives the corrected p-value, corrected across ALL conditions (This is the default column used for all other functions.

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
    covars_renamed  = { c : 'covar_%d' % (i+1) for (i,c) in enumerate(sorted(C.columns, key=len, reverse=True)) }
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
    R['fdr']  = stats.p_adjust(R.pvalue.values, method='fdr')
    return R
#edef

###################################################################################

def summary(diffex, ov=False,
            col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
            col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
            col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh']):
    """
    summary: Returns a summary of the differential expression tests
    Inputs: 
        diffex : The table of tests (e.g. from voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        overlaps: Also return the overlaps with other conditions
            
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output: List of genes significant by the specified conditions
    """
    S = diffex.copy()
    S['significant'] = 1 * (S[col_qval].values < alpha) & (np.abs(S[col_lfc].values) > lfcThresh )
    S['up']          = 1 * ((S.significant == 1) & (S[col_lfc] > 0))
    S['down']        = 1 * ((S.significant == 1) & (S[col_lfc] < 0))
    
    S = S.groupby(col_contr).agg(sum)[['significant','up','down']]
    
    if ov:
        ov = overlaps(diffex, alpha=alpha, lfcThresh=lfcThresh, col_contr=col_contr, col_pval=col_pval, col_qval=col_qval, col_lfc=col_lfc, col_index=col_index)
        S.columns = pd.MultiIndex.from_product([['summary'], S.columns])
        ov.columns = pd.MultiIndex.from_product([['overlaps'], ov.columns])
        return S.join(ov)
    else:
        return S
    #fi
    
#edef

###################################################################################

def overlaps(diffex,
             col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
             col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
             col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
    """
    Overlaps: Returns a table of overlaps the differential expression tests
    Inputs: 
        diffex : The table of tests (e.g. from voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        overlaps: Also return the overlaps with other conditions
            
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output: List of genes significant by the specified conditions
    """
    sigin = significant_in(diffex, alpha=alpha, lfcThresh=lfcThresh, col_contr=col_contr, col_index=col_index,
                           col_pval=col_pval, col_qval=col_qval, col_lfc=col_lfc)
    ov = ops.lst.overlap(list(sigin.values()))
    ov = pd.DataFrame(ov, columns=sigin.keys(), index=sigin.keys())
    return ov
#edef  

###################################################################################

def sigTests(diffex,
             col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
             col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
             col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
    """
    Return the set rows of significant tests.
    Input:
        diffex: the output of a differential expression test (e.g. voom)
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        
        col_qval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        A subset of the significant rows in input diffex
    """
    return diffex[(diffex[col_qval] < alpha) & (diffex[col_lfc].abs() > lfcThresh)]
#edef

###################################################################################

def significant_in(diffex, contr=None, split_updown=False,
                   col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
                   col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
                   col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
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
    
    sigs = sigTests(diffex, alpha=alpha, lfcThresh=lfcThresh, col_contr=col_contr, col_index=col_index,
                    col_pval=col_pval, col_qval=col_qval, col_lfc=col_lfc)
    sigin = {}
    if split_updown:
        sigin = { (c,direction): sigs[(sigs[col_contr] == c) & ((sigs[col_lfc] > 0) if direction == 'up' else (sigs[col_lfc] <= 0) )][col_index].values
                  for c in contr for direction in ['up', 'down']}
    else:
        sigin = { c : sigs[sigs[col_contr] == c][col_index].values for c in contr }
    #fi
    return sigin
#edef

###################################################################################

def gsea_rankings(diffex, rankby='logfc', contr=None, split_updown=False,
                   col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
                   col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
                   col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
    """
    Return a dictionary of significant genes per contrast
    Inputs:
        rankby: 'logfc' | 'pvalue'. Sort by logFC or pvalue
                if logfc  -> all tests with pvalue < alpha are given the lowest ranking
                if pvalue -> all tests 
        diffex: The output of a differential expression (e.g. voom)
        contr: String or list of strings of contrast names (or None)
        split_updown: Split the differential expressions by up/down expressed
    Outputs:
        dict of { contr : (genes, scores) } if split_updown is False
        dict of { contr: { 'up' : (genes, scores), 'down' : (genes, scores) } } if split_updown is True
    """
    if contr is None:
        contr = list(diffex[col_contr].drop_duplicates())
    #fi
    if isinstance(contr, str):
        contr = [ contr ]
    #fi
    
    diffex = diffex.copy()
    diffex = diffex[~diffex[col_index].isna()]
    
    if rankby == 'logfc':
        col_rank = col_lfc
        diffex.loc[diffex[col_qval] >= alpha, col_lfc] = np.nan
        
    elif rankby == 'pvalue':
        col_rank = col_qval
        diffex[col_fdr] = diffex[col_fdr]
        diffex[diffex[col_lfc] >= lfcThresh, col_fdr] = np.nan
    else:
        raise ValueError('Invalid rankby argument.')
    #fi
    
    R = {}
    
    for c in contr:
        DC = diffex[diffex[col_contr] == c]
        if split_updown:
            r_up = -DC[col_rank].values
            r_up[np.isnan(r_up)] = max(r_up[~np.isnan(r_up)])
            
            r_down = DC[col_rank].values
            r_down[np.isnan(r_down)] = max(r_down[~np.isnan(r_down)])
            
            R[(c, 'up')] = (DC[col_index].values, r_up)
            R[(c, 'down')] = (DC[col_index].values, r_down)
        else:
            r = np.abs(DC[col_rank].values)
            r[np.isnan(r)] = max(r[~np.isnan(r)])
            R[c] = (DC[col_index], r)
        #fi
    #efor
    
    return R
#edef

###################################################################################

def isSigIn(diffex, contr=None, logic=all,
            col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
            col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
            col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
    """
    isSigIn: Returns genes that are significant in a given set of contrasts.
    Inputs: 
        diffex : The table of tests
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
        contr = set(diffex[col_contr].values)
    #fi
    relRes = sigTests(diffex, alpha=alpha, lfcThresh=lfcThresh, col_contr=col_contr, col_index=col_index,
                      col_pval=col_pval, col_qval=col_qval, col_lfc=col_lfc)
    sigGeneGroups = ops.lst.group(relRes[[col_contr, col_index]].values, key=lambda x: x[1], value=lambda x: x[0])
    sigGenes = [ g for g in sigGeneGroups if logic([ (c in sigGeneGroups[g]) for c in contr ]) ]
    return sigGenes
#edef

###################################################################################

def volcanoPlot(diffex,contr=None, genes=None, ax=None,
                col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
                col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
                col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'], **kwargs):
    """
    volcanoPlot: Plots volcanoplots for each contrast
    Inputs:
        diffex : The table of tests (e.g. output from voom)
        contr: String or list of Strings. Which contrast to investigate? 
        genes: List[String] | None
            Which genes to plot? None means plot all.
        ax: Matplotlib axes to plot on. Same length as contr
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        
        col_pval: The column to use as pvalue
        col_qval: The column to use as corrected pvalues
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        (Figure, axes)
    """

    tests = diffex[col_contr].drop_duplicates().values

    if contr is not None:
        tests = [ contr ] if isinstance(contr, str) else contr
    #fi
    
    axes = None
    if ax is not None:
        axes = ax if hasattr(ax, '__len__') else [ ax ]
        assert len(ax) == len(tests)
        fig = axes[0].get_figure()
    else:
        ncols = min(np.ceil(np.sqrt(len(tests))), 4)
        nrows = np.ceil(len(tests) / ncols)
        fig, axes = utils.figure.subplots(ncols=int(ncols), nrows=int(nrows),
                                          figsize = (5*ncols, 5*nrows),
                                          sharey=True, sharex=True, **kwargs)
    #fi
    lfc = diffex[col_lfc][np.isfinite(diffex[col_lfc].values)]
    rangeFC   = (lfc.min(), lfc.max())
    rangePV   = -np.log10([diffex[col_pval].max(), diffex[col_pval].min()])
    sigThresh = diffex[diffex[col_qval] < 0.05][col_pval].max()
    sigThresh = -np.log10(sigThresh)
    for idx, test in enumerate(tests):
        relRows = diffex[diffex[col_contr] == test]
        relRows = relRows if genes is None else relRows[relRows[col_index].isin(genes)]
        sigRows = sigTests(relRows, alpha=alpha, lfcThresh=lfcThresh, col_pval=col_pval,
                           col_qval=col_qval, col_lfc=col_lfc, col_contr=col_contr)
        axes[idx].scatter(relRows[col_lfc], -np.log10(relRows[col_pval]), s=0.5, label=None)
        axes[idx].set_title(test)
        axes[idx].plot([rangeFC[0], rangeFC[1]], [sigThresh, sigThresh], linestyle=':', c='red', label='Adjusted pvalue threshold')
        axes[idx].plot([-lfcThresh, -lfcThresh], rangePV, linestyle=':', c='orange')
        axes[idx].plot([lfcThresh, lfcThresh], rangePV, linestyle=':', c='orange', label='logFC threshold')
        axes[idx].legend(loc='upper center')
        axes[idx].scatter(sigRows[col_lfc], -np.log10(sigRows[col_pval]), c='r', s=0.5)
        axes[idx].set_ylabel('-log10 pvalue')
        axes[idx].set_xlabel('log Fold Change')
    #efor
    
    return fig, axes
#edef

###################################################################################

def compair(diffex, contrA, contrB,
            col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
            col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
            col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
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
        col_index: The column to use as the gene identifier
        
    Outputs:
        Dataframe of comparisons, per gene
        
        Comparisons are given:
            significant: Is the gene significantly differentially expressed in any condition (pvalue < alpha) & lfc > lfcThresh
            sig_both:    Is the gene significantly differentially expressed in all conditions?
            consistent:  Does the Log Fold change have the same sign in both contrasts
            up:          Is the LFC consistent and larger than zero in both?
            down:        Is the LFC consistent and less than zero in both?
            stronger:    Is the Log Fold Change consistent and larger in contrast B than in contrast A?
            weaker:      Is the LFC consistent and smaller in contrast B than in contrast A?
            slope:       In terms of the volcano plot, with which slope does the gene move?
            direction:   In terms of LFC, which direction does the gene move? (-1 to left, +1 to right)
    """
    
    D = diffex[diffex[col_contr].isin([contrA, contrB])].copy()
    D = D.pivot(index=col_index, columns=col_contr, values=[col_pval, col_qval, col_lfc])
    D['significant'] = list(map(np.any, (D[col_qval].values < alpha) & (D[col_lfc].abs().values > lfcThresh )))
    D['sig_both']    = list(map(np.all, (D[col_qval].values < alpha) & (D[col_lfc].abs().values > lfcThresh )))
    D['consistent']  = list(map(lambda x: ~np.logical_xor(*x), (D[col_lfc] > 0).values))
    D['up']          = list(map(np.all, (D[col_lfc].values > 0 )))
    D['down']        = list(map(np.all, (D[col_lfc].values < 0 )))
    D['stronger']    = D.consistent.values & (D[col_lfc][contrB].abs().values > D[col_lfc][contrA].abs().values)
    D['weaker']      = D.consistent.values & (D[col_lfc][contrB].abs().values < D[col_lfc][contrA].abs().values)
    D['slope']       = (-np.log10(D[col_pval][contrB].values)+np.log10(D[col_pval][contrA].values))/(D[col_lfc][contrB].values-D[col_lfc][contrA].values)
    D['direction']   = np.sign(D[col_lfc][contrB].values - D[col_lfc][contrA].values)
    
    return D
#edef

###################################################################################

def pairedVolcanoPlot(diffex, contrA, contrB, genes=None, only_significant=True, ax=None, color='stronger',
                      col_contr=_defaults['col_contr'], col_index=_defaults['col_index'],
                      col_pval=_defaults['col_pval'], col_qval=_defaults['col_qval'],
                      col_lfc=_defaults['col_lfc'], alpha=_defaults['alpha'], lfcThresh=_defaults['lfcThresh'] ):
    """
    pairedVolcanoPlot: Plot two volcanoplots on top of each other, and join genes by lines
    Inputs:
        diffex: The table of tests (e.g. output from voom)
        contrA: Contrast A to display (in blue)
        contrB: Contrast B to display (in red)
        genes: List[String] | None
            Which genes to plot. None means plot all.
        only_significant: Only plot genes that are significant in either one of the two contrasts
        color: How to color the plots between the two contrasts?
               See output from compair().
        ax: Matplotlib axes to plot on. Same length as contr
        alpha: pvalue threshold
        lfcThresh: LogFC threshold
        
        
        col_pval: The column to use as corrected pvalue
        col_lfc: The column to use as log fold change
        col_contr: The column to use as the contrast column
    Output:
        List of genes significant by the specified conditions
    """
    
    diffex_sel = diffex if genes is None else diffex[diffex[col_index].isin(genes)]
                                                        
    cmp = compair(diffex_sel, contrA, contrB, alpha=alpha,
                  lfcThresh=lfcThresh, col_contr=col_contr, col_index=col_index,
                  col_pval=col_pval, col_qval=col_qval, col_lfc=col_lfc)
    
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
    
    lfc = diffex[col_lfc][np.isfinite(diffex[col_lfc].values)]
    rangeFC   = (lfc.min(), lfc.max())
    rangePV = -np.log10([diffex[col_pval].min().min(),diffex[col_pval].max().max()])
    sigThresh = diffex[diffex[col_qval] < 0.05][col_pval].max()
    sigThresh = -np.log10(sigThresh)
    
    ax.plot([rangeFC[0], rangeFC[1]], [sigThresh, sigThresh], linestyle=':', c='red', label='Adjusted pvalue threshold', alpha=0.5, linewidth=2)
    ax.plot([-lfcThresh, -lfcThresh], rangePV, linestyle=':', c='orange', alpha=0.5, linewidth=2)
    ax.plot([lfcThresh, lfcThresh], rangePV, linestyle=':', c='orange', alpha=0.5, linewidth=2, label='logFC threshold')

    ax.legend()
    
    ax.set_title('Paired volcano plot\n%s (blue) vs %s (red)' % (contrA, contrB))
    ax.set_xlabel('Log Fold Change')
    ax.set_ylabel('- log10 pvalue')
    
    return ax, cmp
#edef
