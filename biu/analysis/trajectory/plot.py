from ... import utils
from ... import ops

pd        = utils.py.loadExternalModule("pandas")
np        = utils.py.loadExternalModule("numpy")
plt       = utils.py.loadExternalModule("matplotlib.pylab")
Rectangle = utils.py.loadExternalModule("matplotlib.patches", "Rectangle")

from collections import namedtuple

def _trajectoryPlot1D_super(pdata, conditions, component, cmap, ax, pca_var_explained=None, **kwargs):
    """
    Internal function for trajectoryPlot. Plots in 1D
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component: String. which component to plot
        cmap: matplotlib colormap, which colormap to use
        colors: A list of colors for each sample's trajectory. Only works for two conditions...
        pca_var_explained: Dictionary of mapping component -> var_explained
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted=None
    """

    #X = pd.DataFrame(X.values + np.array([x_pos] * X.shape[1]).transpose() + .5, columns=X.columns)
    #Y = pd.DataFrame(Y.values + np.array([y_pos] * Y.shape[1]).transpose() + .5, columns=X.columns)
    X = pdata[component][conditions]
    Y = pd.DataFrame([[ y/2.0 for y in range(X.shape[1]) ] for i in range(X.shape[0])], columns=X.columns)

    for i, cond in enumerate(X.columns):
        ax.scatter(X[cond].values, Y[cond].values, zorder=2, s=1, label=str(cond), c=cmap(i))
    #efor
    ax.legend()
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('%s%s' % (component, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component] * 100)))
    ax.set_ylabel('Conditions')


    lines = ax.plot(X.values.transpose(),
                    Y.values.transpose(),
                    linewidth=0.3, zorder=1, c='k')
    
    return pdata, ax, lines, None
#edef

def _trajectoryPlot1D_grid(pdata, conditions, component, cmap, ax, pca_var_explained=None, reorder=True, reorder_by=None, **kwargs):
    """
    Internal function for trajectoryPlot. Plots in 1D
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component: String. which component to plot
        cmap: matplotlib colormap, which colormap to use
        pca_var_explained: Dictionary of mapping component -> var_explained
        reorder: Boolean, reorder the elements in the grid based on correlation.
        reorder_by: A list of conditions to order by instead of the usual conditions.
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted
    """
    if reorder:
        sorted_pdata = ops.dataframe.reorder(pdata[component][conditions], distance='correlation',
                                             by=(pdata[component][reorder_by] if reorder_by is not None else None))
    else:
        sorted_pdata = pdata[component][conditions]
    #fi
    ylim = (np.min(pdata.values), np.max(pdata.values))
    lines = []
    grid  = []

    for i in range(0, sorted_pdata.shape[0]):
        row = sorted_pdata.iloc[i].values
        xvals = [i]*len(row) + np.arange(0,len(row))/(2*len(row))
        yvals = row
        c = [ cmap(i) for i, w in enumerate(row)]
        ax.scatter(xvals, yvals, zorder=2, c=c)
        l = ax.plot(xvals, yvals, linestyle='-', c='k', zorder=1)
        lines.append(l)
        
        if i%2 == 1:
            rect = Rectangle((np.mean(xvals)-0.5, ylim[0]), 1, ylim[1]-ylim[0],
                             color='gray', alpha=0.2, zorder=1)
            grid.append(rect)
            ax.add_patch(rect)
        #fi
    #efor
    
    custom_lines = [plt.Line2D([0], [0], color=cmap(i), marker='.', lw=0) for i in range(sorted_pdata.shape[1]) ]
    ax.legend(custom_lines, sorted_pdata.columns)
    
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('Parcticipants')
    ax.set_ylabel('%s%s' % (component, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component] * 100)))
    ax.set_ylim(ylim)
    
    return sorted_pdata, ax, lines, grid
#edef

def _trajectoryPlot1D_embed(pdata, conditions, component, cmap, ax, pca_var_explained=None, **kwargs):
    """
    Internal function for trajectoryPlot. Embed 1D measure in 2D-PC space
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component: String. Which component to plot on x-axis
        cmap: matplotlib colormap, which colormap to use
        pca_var_explained: Dictionary of mapping component -> var_explained
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted=None
    """

    max_x = np.max(pdata[component][conditions].values)
    min_x = np.min(pdata[component][conditions].values)


    X = (((pdata[component][conditions] - min_x) / (max_x - min_x))) * 4 - 2
    Y = pd.DataFrame([[ (y/float(X.shape[1])) for y in range(X.shape[1]) ] for i in range(X.shape[0])], columns=X.columns)


    nsamples = pdata.shape[0]
    ncols = int(np.ceil(np.sqrt(nsamples)))
    nrows = int(np.ceil(nsamples / ncols))
    
    embed = ops.dataframe.pca(pdata, nc=2)

    x_pos = embed.PCA_1.values
    y_pos = embed.PCA_2.values

    X = pd.DataFrame(X.values + np.array([x_pos] * X.shape[1]).transpose() + .5, columns=X.columns)
    Y = pd.DataFrame(Y.values + np.array([y_pos] * Y.shape[1]).transpose() + .5, columns=X.columns)


    for i, cond in enumerate(X.columns):
        ax.scatter(X[cond].values, Y[cond].values, zorder=2, s=1, label=str(cond), c=cmap(i))
    #efor
    ax.legend()

    lines = ax.plot(X.values.transpose(),
                    Y.values.transpose(),
                    linewidth=0.3, zorder=1, c='k')
    
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('%s%s' % (component, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component] * 100)))
    
    
    return pdata, ax, lines, None
#edef

def _trajectoryPlot2D_super(pdata, conditions, component1, component2, cmap, ax, pca_var_explained=None, **kwargs):
    """
    Internal function for trajectoryPlot. Plots in 2D
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component1: String. Which component to plot on x-axis
        component2: String. Which component to plot on y-axis
        cmap: matplotlib colormap, which colormap to use
        pca_var_explained: Dictionary of mapping component -> var_explained
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted=None
    """
    for i, cond in enumerate(conditions):
        ax.scatter(pdata[component1][cond],
                   pdata[component2][cond],
                   c=cmap(i), zorder=2, label='%s' % str(cond))

        ax.legend()
    #efor
    
    lines = ax.plot(pdata[component1][conditions].values.transpose(),
                    pdata[component2][conditions].values.transpose(),
                    linewidth=0.3, zorder=1, c='k')    
    
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('%s%s' % (component1, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component1] * 100)))
    ax.set_ylabel('%s%s' % (component2, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component2] * 100)))
    
    return pdata, ax, lines, None
#edef

def _trajectoryPlot2D_grid(pdata, conditions, component1, component2, cmap, ax, pca_var_explained=None, reorder=True, reorder_by=None, **kwargs):
    """
    Internal function for trajectoryPlot. Plot on a grid in 2D
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component1: String. Which component to plot on x-axis
        component2: String. Which component to plot on y-axis
        cmap: matplotlib colormap, which colormap to use
        pca_var_explained: Dictionary of mapping component -> var_explained
        reorder: Order the elements in the grid based on a correlation similarity
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted
    """
    if reorder:
        pdata = ops.dataframe.reorder(pdata, by=(pdata[component][reorder_by if reorder_by is not None else conditions]))
    #fi

    max_x = np.max(pdata[component1][conditions].values)
    min_x = np.min(pdata[component1][conditions].values)
    max_y = np.max(pdata[component2][conditions].values)
    min_y = np.min(pdata[component2][conditions].values)

    X = (((pdata[component1][conditions] - min_x) / (max_x - min_x)) - 0.5)
    Y = (((pdata[component2][conditions] - min_y) / (max_y - min_y)) - 0.5)

    nsamples = pdata.shape[0]
    ncols = int(np.ceil(np.sqrt(nsamples)))
    nrows = int(np.ceil(nsamples / ncols))

    rs = [ Rectangle(xy=(x,y), height=1, width=1,
                     color='gray', alpha=0.2, zorder=1)
           for x in range(ncols) for y in range(nrows) if (x%2==1)^(y%2==1)]
    [ ax.add_patch(r) for r in rs ]
    
    ax.set_xlim([0, ncols])
    ax.set_ylim([0, nrows])

    x_pos = np.array([ x for h in range(nrows) for x in range(ncols) ][:X.shape[0]])
    y_pos = np.array([ y for h in range(nrows) for y in [ h ] * ncols ][::-1][:X.shape[0]])

    X = pd.DataFrame(X.values + np.array([x_pos] * X.shape[1]).transpose() + .5, columns=X.columns)
    Y = pd.DataFrame(Y.values + np.array([y_pos] * Y.shape[1]).transpose() + .5, columns=X.columns)


    for i, cond in enumerate(X.columns):
        ax.scatter(X[cond].values, Y[cond].values, zorder=2, s=1, label=str(cond), c=cmap(i))
    #efor
    ax.legend()

    lines = ax.plot(X.values.transpose(),
                    Y.values.transpose(),
                    linewidth=0.3, zorder=1, c='k')
    
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('%s%s' % (component1, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component1] * 100)))
    ax.set_ylabel('%s%s' % (component2, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component2] * 100)))
    
    return pdata, ax, lines, rs
#edef

def _trajectoryPlot2D_embed(pdata, conditions, component1, component2, cmap, ax, pca_var_explained=None, **kwargs):
    """
    Internal function for trajectoryPlot. Embed in PC space in 2D
    Inputs:
        pdata: DataFrame. A pivoted table with multilevel index on columns: components -> conditions
        conditions: List of strings. conditions to plot
        component1: String. Which component to plot on x-axis
        component2: String. Which component to plot on y-axis
        cmap: matplotlib colormap, which colormap to use
        pca_var_explained: Dictionary of mapping component -> var_explained
        **kwargs: Ignored arguments
    Returns:
        Tuple of plotted_data, axis_plotted_to, lines_plotted, grid_plotted=None
    """

    max_x = np.max(pdata[component1][conditions].values)
    min_x = np.min(pdata[component1][conditions].values)
    max_y = np.max(pdata[component2][conditions].values)
    min_y = np.min(pdata[component2][conditions].values)

    X = (((pdata[component1][conditions] - min_x) / (max_x - min_x)) * 4 - 2)
    Y = (((pdata[component2][conditions] - min_y) / (max_y - min_y)) * 4 - 2)


    nsamples = pdata.shape[0]
    ncols = int(np.ceil(np.sqrt(nsamples)))
    nrows = int(np.ceil(nsamples / ncols))
    
    embed = ops.dataframe.pca(pdata, nc=2)

    x_pos = embed.PCA_1.values
    y_pos = embed.PCA_2.values

    X = pd.DataFrame(X.values + np.array([x_pos] * X.shape[1]).transpose() + .5, columns=X.columns)
    Y = pd.DataFrame(Y.values + np.array([y_pos] * Y.shape[1]).transpose() + .5, columns=X.columns)


    for i, cond in enumerate(X.columns):
        ax.scatter(X[cond].values, Y[cond].values, zorder=2, s=1, label=str(cond), c=cmap(i))
    #efor
    ax.legend()

    lines = ax.plot(X.values.transpose(),
                    Y.values.transpose(),
                    linewidth=0.3, zorder=1, c='k')
    
    ax.set_title('TrajectoryPlot (%s)' % ('->'.join(conditions)))
    ax.set_xlabel('%s%s' % (component1, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component1] * 100)))
    ax.set_ylabel('%s%s' % (component2, '' if pca_var_explained is None else ' (%2.3f)' % (pca_var_explained[component2] * 100)))
    
    
    return pdata, ax, lines, None
#edef

def plot(data, conditions=None, samples=None, dim=2, pca=True,
         component=0, component1=0, component2=1, ax=None, cmap=plt.get_cmap('Set2'),
         display_1d='grid', display_2d='super', reorder=True, reorder_by=None, pca_fit=None):
    """
    Make a 1/2D PC trajectoryplot of samples
    Inputs:
        data: The tdata/pdata object from the Trajectory object
        
        conditions: Which conditions to plot?
                    NOTE: if pca=True, it will perform PCA on ALL conditions, but only plot these.
        samples:    Which samples to plot?
        pca:        Boolean. Do PCA? Otherwise component/component1/component2 refer to columns in the original data matrix
        dim:        Integer [1|2]. How many dimensions to plot.
        col_index:  Which column is the index column?
        col_condition: Which column is the condition column?
        component:  The PC component to use when making a 1D plot
                    if pca=False, then it is either an integer index of the column to use in original data matrix, or the name of the column
        component1: The PC component to use on the x-axis when making a 2D plot
                    if pca=False, then it is either an integer index of the column to use in original data matrix, or the name of the column
        component2: The PC component to use on the y-axis when making a 2D plot
                    if pca=False, then it is either an integer index of the column to use in original data matrix, or the name of the column

        ax: Matplotlib axis object or None. Which axis to plot on?
        cmap: Matplotlib colormap. Which colormap to use
        display_1d: How to display the data in the 2D case.
                    super: Superimpose the plots per index
                    grid:  Plot each index seperately, arranged in a grid
                    embed: Plot each index in an embedded space.
        display_2d: How to display the data in the 2D case.
                    super: Superimpose the plots per index
                    grid:  Plot each index seperately, arranged in a grid
                    embed: Plot each index in an embedded space.
        reorder: Boolean. If possible, order the points based on a hierarchical clustering.
                 Possible in : 2D_grid and 1D_grid
        reorder_by: tdata/pdata/plotted_data from previous trajectory plot. Used to order the elements in the correct way.
                    Must be in same format.
        
    Output:
        A namedtuple, with values:
          * plotted_data: The plotted data
          * ax:           The matplotlib axisn plotted to
          * lines:        The lines plotted
          * grid:         The rectangles plotted on the grid, if relevant, otherwise None
        
    Usage:
        If you have several dataframes, D1, D2, D3... with data on the same samples (in the same order), then:
          trajectoryPlot([D1, D2, D3] ['name_D1', 'name_D2', 'name_D3'])

        If you have one dataframe, D with data, and another dataframe with column information (columns with index and condition), then:
          trajectoryPlot(D, labels)
    """
    
    pca_var_explained = None
    if pca:
        pca_var_explained = { 
          'PCA_%d' % (component+1) : pca_fit.explained_variance_ratio_[component],
          'PCA_%d' % (component1+1) : pca_fit.explained_variance_ratio_[component1],
          'PCA_%d' % (component2+1) : pca_fit.explained_variance_ratio_[component2]}
        component  = 'PCA_%d' % (component+1)
        component1 = 'PCA_%d' % (component1+1)
        component2 = 'PCA_%d' % (component2+1)
    else:
        component  = data.columns[component] if (component not in data.columns) and (isinstance(component, int)) else component
        component1 = data.columns[component1] if (component1 not in data.columns) and (isinstance(component1, int)) else component1
        component2 = data.columns[component2] if (component2 not in data.columns) and (isinstance(component2, int)) else component2
    #fi

    if ax is None:
        fig, axes = utils.figure.subplots(ncols=1, nrows=1)
        ax = axes[0]
    #fi
    
    if dim == 1:
        funcs = { 'super' : _trajectoryPlot1D_super,
                  'grid' :  _trajectoryPlot1D_grid,
                  'embed':  _trajectoryPlot1D_embed}
        p_data, ax, lines, grid = funcs[display_1d.lower()](pdata=data, conditions=conditions, component=component,
                                                            cmap=cmap, ax=ax, pca_var_explained=pca_var_explained,
                                                            reorder=reorder, reorder_by=reorder_by)
    else:
        funcs = { 'super' : _trajectoryPlot2D_super,
                  'grid':   _trajectoryPlot2D_grid,
                  'embed' : _trajectoryPlot2D_embed}
        p_data, ax, lines, grid = funcs[display_2d.lower()](pdata=data, conditions=conditions,
                                                            component1=component1, component2=component2,
                                                            cmap=cmap, ax=ax, pca_var_explained=pca_var_explained,
                                                            reorder=reorder, reorder_by=reorder_by)
    #fi
    
    nt = namedtuple("Trajectory_plot_data", [ 'pdata', 'ax', 'lines', 'grid' ])
    return nt(p_data, ax, lines, grid)
#edef
         