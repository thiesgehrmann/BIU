from .. import utils

def manhattanPlot(er, chrs=None, ax=None, text=True,
                  col_chr='chromosome', col_pos='position', col_pval='frequentist_add_pvalue',
                  col_annot='phenotype', annotate=False, annots=None, threshold=1e-8):
    """
    
    """
    if chrs is None:
        chrs = sorted(list(set(er[col_chr].values)))
    #fi

    if ax is None:
        fig, axes = utils.figure.subplots(dpi=1000, figsize=(6,2))
        ax = axes[0]
    #fi
                  
    colors = ['#a6611a','#dfc27d','#80cdc1','#018571',]
    ncolors = len(colors)
    er['logp'] = er[col_pval].apply(lambda x: -np.log10(x))
    max_lpval = max(er['logp'])
    min_lpval = min(er['logp'])
    lastpos = 0
    for i, chrID in enumerate(chrs):
        relE = er[er[col_chr] == chrID]
        ax.scatter(lastpos + relE[col_pos], relE['logp'], alpha=0.5, s=0.5, c=colors[i % ncolors])
        ax.set_ylabel('-log10(pvalue)')
        nextpos = lastpos + max(relE[col_pos])
        if text:
            ax.text( (nextpos+lastpos)/2, min_lpval-1, str(chrID), horizontalalignment='center', fontsize=6)
        #fi
        
        if annotate or (annots is not None):
            sigE = relE[(relE[col_pval] <= threshold)]
            if annots is not None:
                sigE = sigE[(sigE.phenotype.isin(annots))]
            #fi
            
            grouped = sigE[sigE.groupby(col_annot).transform(max).logp == sigE.logp].drop_duplicates(col_annot)
            for i, gene in grouped.iterrows():
                ax.text(lastpos+gene['position'], gene['logp'], gene[col_annot], fontsize=4)
        #fi
        
        lastpos = nextpos
    #efor
    ax.get_xaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xlim((0, lastpos))
    ax.plot([0, lastpos], -np.log10([threshold,threshold]), c='r', linewidth=0.5)
    
    return ax
#edef

def miamiPlot(er1, er2, **kwargs):
    fig, axes = utils.figure.subplots(dpi=1000, figsize=(6,4), ncols=1, nrows=2)
    manhattanPlot(er1, ax=axes[0], text=True, **kwargs)
    manhattanPlot(er2, ax=axes[1], text=False, **kwargs)
    
    min_ylim = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
    max_ylim = max(axes[1].get_ylim()[0], axes[1].get_ylim()[1])
    axes[0].set_ylim(min_ylim, max_ylim)
    axes[1].set_ylim(min_ylim, max_ylim)
    
    axes[1].invert_yaxis()
    
    return fig, axes
#edef
