from .. import pyUtils as py

import math
mpl = py.loadExternalModule('matplotlib')
plt = py.loadExternalModule('matplotlib', 'pylab')
sns = py.loadExternalModule('seaborn')
sstats = py.loadExternalModule('scipy', 'stats')
pd  = py.loadExternalModule('pandas')
np  = py.loadExternalModule('numpy')

__all__ = [ 'SingleJointPlot', 'MultiJointPlot', 'jointplot', 'multijointplot' ]

class SingleJointPlot(object):
    """
    The grid for a SingleJointPlot.
    
    You can generate a SingleJointPlot by:
    
    sjp = SingleJointPlot()
    sjp.plot(x, y)
    
    """
    def __init__(self, hist_x=None, hist_y=None, joint=None):
        """Generate a singleJointPlot grid
        
        parameters:
        -----------
        hist_x, hist_y, joint: matplotlib axes objects
            The axes to which we will plot
            If they are all None, then we will create them automatically.
    
        Output:
        SingleJointPlot object
        
        """
        
        if (hist_x is None) and (hist_y is None) and (joint is None):
            mjp = MultiJointPlot(1)
            hist_x, hist_y, joint = mjp.axes[0]
        #fi
        
        self.hist_x = hist_x
        self.hist_y = hist_y
        self.joint  = joint
    #edef
    
    
    def plot(self, x, y, data=None, kind='scatter', dropna=True,
             title=None, xlabel=None, ylabel=None, color=None,
             joint_kws=None, marginal_kws=None, annot_kws=None):
        """Draw a plot of two variables with bivariate and univariate graphs.

        parameters:
        ----------
        x, y : strings or vectors
            Data or names of variables in ``data``.
        data : DataFrame, optional
            DataFrame when ``x`` and ``y`` are variable names.
        kind : { "scatter" | "reg" | "resid" | "kde" | "hex" }, optional
            Kind of plot to draw.
        color : matplotlib color, optional
            Color used for the plot elements.
        dropna : bool, optional
            If True, remove observations that are missing from ``x`` and ``y``.
        title: String, optional
            String to set as title
        xlabel: String, optional
            String to set as x label on joint axis
        ylabel: String, optional
            String to set as y label on joint axis
        {joint, marginal, annot}_kws : dicts, optional
            Additional keyword arguments for the plot components.
            
        Output
        -------
        The SingleJointPlot which was plotted to
        """
        
        joint_kws = {} if joint_kws is None else joint_kws
        marginal_kws = {} if marginal_kws is None else marginal_kws
        annot_kws = {} if joint_kws is None else annot_kws

        if data is not None:
            xlabel = x if xlabel is None else xlabel
            ylabel = y if ylabel is None else ylabel
            title  = '%s / %s' % (x, y) if title is None else title
            x = data[x].values
            y = data[y].values
        else:
            x = np.array(x)
            y = np.array(y)
            xlabel = 'X' if xlabel is None else xlabel
            ylabel = 'Y' if ylabel is None else ylabel
            title  = '' if title is None else title
        #fi
        
        if dropna:
            keep = ~pd.isna(x) & ~pd.isna(y) 
            x = x[keep]
            y = y[keep]
        #fi
        
        if color is None:
            color = sns.palettes.color_palette()[0]
        #fi
        color_rgb = mpl.colors.colorConverter.to_rgb(color)
        colors = [sns.utils.set_hls_values(color_rgb, l=l)  # noqa
                  for l in np.linspace(1, 0, 12)]
        cmap = sns.palettes.blend_palette(colors, as_cmap=True)

        if kind == "scatter":
            joint_kws.setdefault("color", color)
            self._plot_joint(x, y, plt.scatter, **joint_kws)

            marginal_kws.setdefault("kde", False)
            marginal_kws.setdefault("color", color)
            self._plot_marginals(x, y, sns.distplot, **marginal_kws)

        elif kind.startswith("hex"):

            x_bins = min(sns.distributions._freedman_diaconis_bins(x), 50)
            y_bins = min(sns.distributions._freedman_diaconis_bins(y), 50)
            gridsize = int(np.mean([x_bins, y_bins]))

            joint_kws.setdefault("gridsize", gridsize)
            joint_kws.setdefault("cmap", cmap)
            self._plot_joint(x, y, plt.hexbin, **joint_kws)

            marginal_kws.setdefault("kde", False)
            marginal_kws.setdefault("color", color)
            self._plot_marginals(x, y, sns.distplot, **marginal_kws)

        elif kind.startswith("kde"):

            joint_kws.setdefault("shade", True)
            joint_kws.setdefault("cmap", cmap)
            self._plot_joint(x, y, sns.kdeplot, **joint_kws)

            marginal_kws.setdefault("shade", True)
            marginal_kws.setdefault("color", color)
            self._plot_marginals(x, y, sns.kdeplot, **marginal_kws)

        elif kind.startswith("reg"):

            marginal_kws.setdefault("color", color)
            self._plot_marginals(x, y, sns.distplot, **marginal_kws)

            joint_kws.setdefault("color", color)
            self._plot_joint(x, y, sns.regplot, **joint_kws)

        elif kind.startswith("resid"):

            joint_kws.setdefault("color", color)
            
            self._plot_joint(x, y, sns.residplot, **joint_kws)

            marginal_kws.setdefault("color", color)
            marginal_kws.setdefault("kde", False)
            sns.distplot(x, ax=self.hist_x, **marginal_kws)
            sns.distplot(y, vertical=True, fit=sstats.norm, ax=self.hist_y,
                     **marginal_kws)
            stat_func = None
        else:
            msg = "kind must be either 'scatter', 'reg', 'resid', 'kde', or 'hex'"
            raise ValueError(msg)
        #fi

        (self.joint if (self.hist_x is None) else self.hist_x).set_title(title)
        self.joint.set_xlabel(xlabel)
        self.joint.set_ylabel(ylabel)
        
        return self
    #edef
    
    def _plot_joint(self, x, y, func, **kwargs):
        """Draw a bivariate plot of `x` and `y`.
        Parameters
        ----------
        func : plotting callable
            This must take two 1d arrays of data as the first two
            positional arguments, and it must plot on the "current" axes.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.
        Returns
        -------
        self : JointGrid instance
            Returns `self`.
        """
        plt.sca(self.joint)
        func(x, y, **kwargs)

        return self

    def _plot_marginals(self, x, y, func, **kwargs):
        """Draw univariate plots for `x` and `y` separately.
        Parameters
        ----------
        func : plotting callable
            This must take a 1d array of data as the first positional
            argument, it must plot on the "current" axes, and it must
            accept a "vertical" keyword argument to orient the measure
            dimension of the plot vertically.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.
        Returns
        -------
        self : JointGrid instance
            Returns `self`.
        """
        kwargs["vertical"] = False
        plt.sca(self.hist_x)
        func(x, **kwargs)

        kwargs["vertical"] = True
        plt.sca(self.hist_y)
        func(y, **kwargs)
    #edef

#eclass

class MultiJointPlot(object):
    """
    Make a grid for a multijointplot.
    
    Has some objects:
    * n_pairs: (Integer) Number of plots possible
    * axes: (List of SingleJointPlots) axes for plotting.
    * plot_axes: (List of matplotlib axes) The axes of the center plots
    * hist_x_axes: (List of matplotlib axes) The axes of the x histograms
    * hist_y_axes: (List of matplotlib axes) The axes of the y histograms
    * _dims: (Dict) Details about the structure of the grid
    
    Usage:
    
    mjp = MultiJointPlot(3, nrows=1)
    
    mjp.axes[0].plot(x, y)
    mjp.axes[1].plot(a, b, kind='hex')
    mjp.axes[2].plot(c, d, kind='scatter')
    
    
    """
    def __init__(self, n_pairs, nrows=None, ncols=None,
                 hist_part=2, gap_part=0, plot_part=7, gap_between_plots=1,
                 xlim=None, ylim=None,
                 **kwargs):
        """
        Make a grid for a multijointplot
        Parameters:
        -----------
        
        n_pairs: Integer
            How many plots should there by
        nrows, ncols: Integers, optional
            How many rows and columns should there be?
        hist_part: Integer, optional
            Per plot, how many parts of vertical height should be dedicated to the histogram?
        gap_part: Integer, optional
            Per plot, how many parts of vertical height should be dedicated to the gap between histogram and plot?
        plot_part: Integer, optional
            Per plot, how many parts of vertical height should be dedicated to the plot?
        gap_between_plots: Integer
            How many parts of vertical/horizontal space should be between plots?
        {x,y}lim: 2-tuple of floats
            Pre-specified limits to the axes X and Y axes.
            
        **kwargs: Dict
            Additional arguments to matplotlib.figure()
            
        Output:
        A MultiJointPlot object
        """
        
        self.__n_pairs = n_pairs
        
        square = hist_part + gap_part + plot_part
        
        if (nrows is None) and (ncols is not None):
            nrows = math.ceil(n_pairs / ncols)
        elif (ncols is None)  and (nrows is not None):
            ncols = math.ceil(n_pairs / nrows)
        elif (ncols is None) and (nrows is None):
            nrows = math.floor(math.sqrt(n_pairs))
            ncols = math.ceil(n_pairs / nrows)
        #fi
        
        total_grid = (square * nrows + gap_between_plots*(nrows-1),
                      square * ncols + gap_between_plots*(ncols-1))

        fig = plt.figure(figsize=(total_grid[1], total_grid[0]))
        gs =  mpl.gridspec.GridSpec(total_grid[0], total_grid[1])
        
        hist_x_pos = (0,0)
        hist_x_dim = (square - (hist_part + gap_part), hist_part)
        
        hist_y_pos = (plot_part+gap_part, hist_part+gap_part)
        hist_y_dim = (hist_part, square - (hist_part + gap_part))
        
        plot_pos = (0, hist_part + gap_part)
        plot_dim = (plot_part, plot_part)
        
        self._dims = { 'hist_part' : hist_part,
                       'gap_part'  : gap_part,
                       'plot_part' : plot_part,
                       'gap_between_squares' : gap_between_plots,
                       'square'    : square,
                       'nrows' : nrows,
                       'ncols' : ncols,
                       'total_grid' : total_grid}

        self.figure = fig

        self.axes = []
        pair = 0
        for row in range(nrows):
            for col in range(ncols):
                
                if pair >= self.__n_pairs:
                    break
                #fi
                
                base_x = col * (square + gap_between_plots)
                base_y = row * (square + gap_between_plots)

                hist_x = None
                hist_y = None
                if hist_part > 0:
                    hist_x = fig.add_subplot(gs[ base_y+hist_x_pos[1]:base_y+hist_x_pos[1]+hist_x_dim[1],
                                                 base_x+hist_x_pos[0]:base_x+hist_x_pos[0]+hist_x_dim[0] ])
                    
                    # Turn off tick labels
                    hist_x.set_xticklabels([])
                    hist_x.spines['right'].set_visible(False)
                    hist_x.spines['top'].set_visible(False)
                    hist_x.spines['bottom'].set_visible(False)
                    if xlim is not None:
                        hist_x.set_xlim(xlim)
                    #fi

                    hist_y = fig.add_subplot(gs[ base_y+hist_y_pos[1]:base_y+hist_y_pos[1]+hist_y_dim[1],
                                                 base_x+hist_y_pos[0]:base_x+hist_y_pos[0]+hist_y_dim[0] ])
                    
                    hist_y.set_yticklabels([])
                    hist_y.spines['right'].set_visible(False)
                    hist_y.spines['top'].set_visible(False)
                    hist_y.spines['left'].set_visible(False)
                    if ylim is not None:
                        hist_y.set_xlim(ylim)
                    #fi
                else:
                    hist_x = hist_y = None
                #fi

                joint = fig.add_subplot(gs[ base_y+plot_pos[1]:base_y+plot_pos[1]+plot_dim[1],
                                           base_x+plot_pos[0]:base_x+plot_pos[0]+plot_dim[0] ])
                if xlim is not None:
                    joint.set_xlim(xlim)
                #fi
                if ylim is not None:
                    joint.set_ylim(ylim)
                #fi

                self.axes.append(SingleJointPlot(hist_x, hist_y, joint))
                
                pair = pair + 1
            #efor
        #efor
    #edef
    
    @property
    def joint_axes(self):
        """All center plot axes"""
        return [ sjp.joint for sjp in self.axes ]
    #edef
    
    @property
    def hist_x_axes(self):
        """All X histogram plot axes"""
        return [ sjp.hist_x for sjp in self.axes ]
    #edef
    
    @property
    def hist_y_axes(self):
        """All Y histogram plot axes"""
        return [ sjp.hist_y for sjp in self.axes ]
    #edef
    
#eclass

def jointplot(*pargs, **kwargs):
    """
    A wrapper for `SingleJointPlot.plot()` See the docstring for this function.
    
    Returns:
    MultiJointPlot object with a single plot
    """
    mjp = MultiJointPlot(1)
    
    mjp.axes[0].plot(*pargs, **kwargs)
    
    return mjp
#edef

def multijointplot(pairs=None, data=None, variables=None, kind='scatter', dropna=True, color=None,
                   joint_kws={}, marginal_kws={}, annot_kws={}, **kwargs):
    """
    Automatically plot several pairs of data with jointplots
    
    Parameters:
    pairs: List of 2-tuples OR dictionary with 2-tuple values, optional
        The 2-tuples may either be pairs of numeric arrays to be plotted (x,y), OR
        They may be identifiers referring to column names in a dataframe `data`
    data: pandas dataframe, optional
        If `data` is specified, it is assumed that whatever is provided in `pairs` refers to dataframe columns.
        If `pairs` is unspecified, then all pairs of dataframe columns are generated
    variables: List of pandas dataframe column identifiers
        If `pairs` is unspecified, then all pairs will be generated from the variable names in this list
    kind : { "scatter" | "reg" | "resid" | "kde" | "hex" }, optional
            Kind of plot to draw.
    color : matplotlib color, optional
        Color used for the plot elements.
    dropna : bool, optional
        If True, remove observations that are missing from ``x`` and ``y``.
    {joint, marginal, annot}_kws : dicts, optional
        Additional keyword arguments for the plot components.
    **kwargs: additional arguments, optional
        Additional arguments given to MultiJointPlot
        
    Output:
    A MultiJointPlot Object with plotted figures.    
    """

    if data is not None:
        variables = data.columns if variables is None else variables
        if pairs is None:
            pairs = [ (a,b) for a in variables for b in variables if a != b ]
        #fi
        
        n_pairs = len(pairs)
        
        V_data = { (v1,v2) : (data[v1].values,data[v2].values) for (v1,v2) in pairs }
    elif pairs is not None:
        if isinstance(data, dict):
            V_data = data
        else:
            V_data = { ('A', 'B') % (i+1) : (x,y) for i, (x,y) in enumerate(data) }
        #fi
    else:
        raise ValueError("You must specify x and y as arrays of data, or as strings, together with a dataframe (data)")
    #fi
    
    mjp = MultiJointPlot(n_pairs, **kwargs)
    
    for pair_i, pair in enumerate(V_data):
        x, y = V_data[pair]
        mjp.axes[pair_i].plot(x, y, 
                  title=' / '.join(pair),
                  xlabel=pair[0], ylabel=pair[1],
                  color=color)
    #efor
    
    return mjp
#edef        