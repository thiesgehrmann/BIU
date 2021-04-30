from ... import utils
from ... import ops

np  = utils.py.loadExternalModule('numpy')
plt = utils.py.loadExternalModule('matplotlib.pylab')

###############################################################################

def stacks(D, ax=None, cmap='tab20', reorder=True, legend=True, plot_min_value=False):
    """
    stacks: Plot a stacked barchart from a dataframe.
    
    parameters:
    ------------
    
    D: pandas.Dataframe
        Columns are taxa
        Rows are samples
        You should make sure that rows adds up to 100% YOURSELF
    ax: matplotlib.pylab.Axis
        Axis to plot on. if None, new axis will be made
    cmap: string|callable
        Colormap to use
    reorder: bool
        Reorder the rows based on a bray-curtis complete hierarchical clustering
    legend: bool
        Plot the legend (column names)
    plot_min_value: bool
        This is useful if you have a pseudocount.
    """
    if ax is None:
        fig, axes = utils.figure.subplots(figsize=(15,5))
        ax = axes[0]
    #fi
    
    if not hasattr(cmap, '__call__'):
        cmap = plt.get_cmap(cmap)
    #fi
    
    if reorder:
        D = ops.dataframe.reorder(D, distance='braycurtis', method='complete')
    #fi
    
    mv = D.min().min()
    
    for i, (index, row) in enumerate(D.iterrows()):
        print('\r%d/%d' %(i+1, len(D)), end='')
        cumsum = np.cumsum(row)
        for j, (height, ypos) in enumerate(zip(row, cumsum)):
            if (height == mv) & (not plot_min_value):
                continue
            #fi
            rectangle = plt.Rectangle((i,ypos-height), 1, height, fc=cmap(j), ec=None)
            ax.add_patch(rectangle)
        #efor
    #efor
    
    if legend:
        legend_elements = [plt.Rectangle((0,0), 1, 1, fc=cmap(i), ec=None, label=c) for (i,c) in enumerate(D.columns) ]
        ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    #fi
    
    ax.set_xlim(0,D.shape[0])
    ax.set_ylim(0,1)
    return ax
#edef