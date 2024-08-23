from ... import utils
from ... import ops

np  = utils.py.loadExternalModule('numpy')
plt = utils.py.loadExternalModule('matplotlib.pylab')

###############################################################################

def stacks(D, ax=None, cmap='tab20', hatch=None, reorder=True, legend=True, plot_min_value=False):
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
    
    import shapely.geometry as sg
    import shapely.ops as so
    from matplotlib.patches import Polygon
    import matplotlib as mpl
    
    if ax is None:
        fig, axes = utils.figure.subplots(figsize=(15,5))
        ax = axes[0]
    #fi
    
    if not hasattr(cmap, '__call__'):
        cmap = plt.get_cmap(cmap)
    #fi
    if hatch is None:
        hatch = lambda x: None
    elif not hasattr(hatch, '__call__'):
        hatch_input = hatch
        hatch = lambda x: hatch_input[x]
    #fi
    
    if reorder:
        D = ops.dataframe.reorder(D, distance='braycurtis', method='complete')
    #fi
    
    mv = D.min(axis=1).max()
    
    patches = { j : [ sg.box(0, -0.0001, D.shape[0], 0) ] for j,c in enumerate(D.columns) }
    
    for i, (index, row) in enumerate(D.iterrows()):
        print('\r%d/%d' %(i+1, len(D)), end='')
        cumsum = np.cumsum(row)
        for j, (height, ypos) in enumerate(zip(row, cumsum)):
            if (height == mv) & (not plot_min_value):
                continue
            #fi
            #rectangle = plt.Rectangle((i,ypos-height), 1, height, fc=cmap(j), ec=None)
            box = sg.box(i-0.0001, -0.0001,i+1, ypos)
            patches[j] = patches.get(j,[]) + [box]
            #ax.add_patch(rectangle)
        #efor
    #efor

    
    def select_hatch_color(background_color, light_color="#FFFFFF", dark_color="#000000"):
        rgb = [ c / 255 for c in mpl.colors.to_rgb(background_color) ]
        rgb = [ c / 12.92 if c <= 0.03928 else np.power((c + 0.055) / 1.055, 2.4) for c in rgb ]
        L = (0.2126 * rgb[0]) + (0.7152 * rgb[1]) + (0.0722 * rgb[2]);
        return  dark_color if (L > 0.179) else light_color
    #edef
    
    def select_hatch_color(background_color, light_color="#FFFFFF", dark_color="#000000"):
        r,g,b = [ c * 255 for c in mpl.colors.to_rgb(background_color) ]
        score = (r*0.299 + g*0.587 + b*0.114)
        if score > 186:
            return "#000000"
        else:
            return "#ffffff"
        #fi
    #edef
    
    
    legend_elements = []
    
    for j, patch_group in enumerate(patches):
        merged_patch = so.cascaded_union(patches[patch_group])
        xs, ys = merged_patch.exterior.xy
        mpl.rcParams['hatch.linewidth'] = 0.5
        mpl.rcParams['hatch.color'] = select_hatch_color(cmap(j))
        p = Polygon(list(zip(xs,ys)), fc=cmap(j), ec=None, hatch=hatch(j), zorder=len(patches)-j)
        ax.add_patch(p)
        #break
        #ax.fill(xs, ys, fc=cmap(j), ec=None, zorder=len(patches)-j)
        
        if legend:
            legend_elements.append(plt.Rectangle((0,0), 1, 1, fc=cmap(j), hatch=hatch(j), ec=None, label=D.columns[j]))
        #fi
    #efor
    
    if legend:
        ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4)
    #fi
    
    ax.set_xlim(0,D.shape[0])
    ax.set_ylim(0,1)
    return ax
#edef