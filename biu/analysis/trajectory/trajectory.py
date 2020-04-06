from ... import utils
from ... import ops

pd        = utils.py.loadExternalModule("pandas")
np        = utils.py.loadExternalModule("numpy")
plt       = utils.py.loadExternalModule("matplotlib.pylab")

from .plot import plot as traj_plot

from collections import namedtuple

class Trajectory(object):
    """
    A class to perform trajectory analysis
    """
    def __init__(self, data, labels=None, col_index='index', col_condition='condition', max_components=5):
        """
        Prepare data for a trajectory analysis
        Inputs:
            data, labels: Either-
              * data:   A list of dataframes where columns are measurements and rows are samples
                        Rows and columns are assumed to be the same order in all dataframes
                        Each dataframe is a different condition
                labels: List of Labels for the points (or None, in which case they are named p1, p2, ..., pn)

              * data:   The dataframe of expression/measurements. Rows are samples, Genes/Measurements are columns
                labels: The dataframe of information about the samples. Must have same index as data dataframe.
                        Rows are samples, columns are measurements.
                col_index: String. The index to use in labels for the sample groupings (e.g. individuals)
                col_condition: String. The index to use in labels for condition grouping.
        Outputs:
            Trajectory Object
        """

        if not(isinstance(data, pd.DataFrame) and isinstance(labels, pd.DataFrame)):
            if labels is None:
                labels = [ 'p%d' % (i+1) for i in range(len(data)) ]
            #fi
            conditions = labels
            labels = pd.DataFrame.from_dict({"condition":[ labels[i] for i, p in enumerate(data) for _ in range(p.shape[0])],
                                             "index": list(range(data[0].shape[0])) * len(data)})
            data   = pd.concat([p.reset_index(drop=True) for p in data]).reset_index(drop=True)
            col_condition = 'condition'
            col_index     = 'index'
        #fi
        
        data = data.loc[labels.index]

        pca_data, pca_fit = ops.dataframe.pca(data, nc=max_components, ret_fit=True)

        self.labels = labels.copy()
        self.rdata = data.copy()
        self.tdata = pd.pivot_table(data.join(labels[[col_index, col_condition]]), index=col_index, columns=col_condition, values=list(data.columns))
        self.pdata = pd.pivot_table(pca_data.join(labels[[col_index, col_condition]]), index=col_index, columns=col_condition, values=list(pca_data.columns))
        self.pca_fit = pca_fit
        
        self._meas       = list(self.rdata.columns)
        self._samples    = list(self.tdata.index)
        self._conditions = list(sorted(set(self.labels[col_condition].values)))
    #edef
    
    def plot(self, conditions=None, samples=None, dim=2, pca=True,
             component=0, component1=0, component2=1, ax=None, cmap=plt.get_cmap('Set2'),
             display_1d='grid', display_2d='super', reorder=True, reorder_by=None, pca_fit=None):
        """
        Make a 1/2D PC trajectoryplot of samples
        Inputs:
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
              * plotted_data:   The plotted data
              * ax:      The matplotlib axisn plotted to

        Usage:
            If you have several dataframes, D1, D2, D3... with data on the same samples (in the same order), then:
              trajectoryPlot([D1, D2, D3] ['name_D1', 'name_D2', 'name_D3'])

            If you have one dataframe, D with data, and another dataframe with column information (columns with index and condition), then:
              trajectoryPlot(D, labels)
        """
        
        if conditions is None:
            conditions = self._conditions
        #fi

        if samples is None:
            samples = self._samples
        #fi

        return traj_plot(self.pdata if pca else self.tdata, conditions=conditions, samples=samples, dim=dim,
                         pca=pca, component=component, component1=component1, component2=component2,
                         ax=ax, cmap=cmap, display_1d=display_1d, display_2d=display_2d,
                         reorder=reorder, reorder_by=reorder_by, pca_fit=self.pca_fit)
    #edef
    
    @utils.decorators.deprecated("axis_of_effect has been replaced with mean_axis_of_effect_2.")
    def axis_of_effect(self, cond1, cond2, components=None):
        """
        Determine the trajectory of the mean axis of effect between two conditions, per individual
        Intputs:
            cond1: The name of condition 1
            cond2: The name of condition 2
            components: List of Integers, or None.
                        If you only want to investigate the trajectory in a specific set of components,
                        specify them here. Starts at 0
        Outputs:
            Pandas Series of axis vector, per individual
        """
        if components is None:
            components = list(range(len(self.pdata.columns.levels[0])))
        #fi
        components     = np.array(components, dtype='int')
        not_components = np.array([ c for c in range(len(self.pdata.columns.levels[0])) if c not in components ], dtype='int') 

        pos_1 = self.pdata.swaplevel(0,1, axis=1)[cond1].values
        pos_2 = self.pdata.swaplevel(0,1, axis=1)[cond2].values

        vector = pos_2 - pos_1
        vector[:, not_components] = 0
        return pd.DataFrame(vector, index=self.pdata.index, columns=self.pdata.columns.levels[0])
    #edef
        
    class Mean_Axis_Of_Effect(object):
        """
        Contains the mean axis of effect in a trajectory analysis.
        """
        def __init__(self, center_1, center_2, components, cond_1, cond_2):
            """
            Initialize a Mean_Axis_Of_Effect object.

            parameters:
            -----------
            center_1: numpy array. Mean PC of condition 1
            center_2: numpy array. Mean PC of condition 2
            components: numpy array: Which components of the PCA to use.
            cond_1: String. Name of the first condition
            cond_2: String. Name of the second condition

            All components not listed in `components` will be set to zero.
            """

            self._components = components
            self._not_components = [ c for c in range(len(center_1)) if c not in self._components ]

            self._center_1 = center_1
            self._center_2 = center_2

            self._center_1[self._not_components] = 0
            self._center_2[self._not_components] = 0

            self._vector = self._center_2 - self._center_1

            self._conditions = (cond_1, cond_2)
        #edef

        @property
        def center_1(self):
            """
            Return the first center.
            """
            return self._center_1
        #edef

        @property
        def center_2(self):
            """
            Return the second center.
            """
            return self._center_2
        #edef

        @property
        def vector(self):
            """
            Return the difference between the two centers
            """
            return self._vector
        #edef

        @property
        def conditions(self):
            """
            Return a tuple of the names of the first and second conditions.
            """
            return self.conditions
        #edef

        def plot(self, ax, x_dim=0, y_dim=1, scatter_label='Center', plot_label='Axis of effect', **kwargs):
            """
            Plot the mean axis in two dimensions to an axis.

            parameters:
            -----------
            ax: Matplotlib axis. The axis to plot onto.
            c:  Color. The color of the mean points and lines
            x_dim: Which component to use as the X dimension (of the provided components)
            y_dim: Which component to use as the Y dimension (of the provided components)
            scatter_label: A label for the centers
            plot_label: A label for the axis of effect line.

            **kwargs: Additional arguments to ax.scatter and ax.plot (e.g. color, alpha)

            Returns:
            --------
            matplotlib axis.
            """
            x = self._components[x_dim]
            y = self._components[y_dim]

            c1 = self.center_1
            c2 = self.center_2

            slope_2d  = self.vector[y] / self.vector[x]
            intercept_2d = ((c1[y] + c2[y]) - slope_2d*(c1[x] + c2[x])) / 2


            xrange = ax.get_xlim()

            ax.scatter([c1[x], c2[x]], [c1[y], c2[y]], zorder=2, label=scatter_label, **kwargs)
            ax.plot(xrange, slope_2d * np.array(xrange) + intercept_2d, zorder=2, label=plot_label, **kwargs)
            ax.legend()

            return ax
        #edef

        def __str__(self):
            """
            A string representation of the object
            """
            dstr = "Mean Axis Of Effect object\n\n"
            dstr += "Axis between condition '%s' and '%s'\n" % self._conditions
            dstr += "Components selected: %s\n" % ', '.join(self._components)
            return dstr
        #edef

        def __repr__(self):
            return str(self)
        #edef
    #eclass
    
    def mean_axis_of_effect2(self, cond1, cond2, components=None, ax=None, samples=None, **kwargs):
        """
        Determine the trajectory of the mean axis of effect between two conditions
        Intputs:
            cond1: The name of condition 1
            cond2: The name of condition 2
            components: List of Integers, or None.
                        If you only want to investigate the trajectory in a specific set of components,
                        specify them here. Starts at 0
            ax: Matplotlib axis. if not None, then plot the axis here.
            samples: The samples you want to use
            c:  Matplotlib color. If plotting, then draw with this color
        Outputs:
            Mean_Axis_Of_Effect object
        """
        
        if components is None:
            components = list(range(len(self.pdata.columns.levels[0])))
        #fi
        components = np.array(components, dtype='int')
        
        data = self.pdata if samples is None else self.pdata.loc[samples]

        center_1 = data.swaplevel(0,1, axis=1).mean()[cond1].values
        center_2 = data.swaplevel(0,1, axis=1).mean()[cond2].values
        
        return self.Mean_Axis_Of_Effect(center_1, center_2, components, cond1, cond2)
    #edef
        

    @utils.decorators.deprecated("mean_axis_of_effect has been replaced with mean_axis_of_effect_2.")
    def mean_axis_of_effect(self, cond1, cond2, components=None, ax=None, c='b'):
        """
        Determine the trajectory of the mean axis of effect between two conditions
        Intputs:
            cond1: The name of condition 1
            cond2: The name of condition 2
            components: List of Integers, or None.
                        If you only want to investigate the trajectory in a specific set of components,
                        specify them here. Starts at 0
            ax: Matplotlib axis. if not None, then plot the axis here.
            c:  Matplotlib color. If plotting, then draw with this color
        Outputs:
            namedtuple with:
                slope_2d:     in the first two components, the slope
                intercept_2d: In the first two components, the intercept
                mean_1_2d:    in the first two components, the mean of condition1
                mean_2_2d:    in the first two components, the mean of condition2
                vector:       The difference between the mean of cond1 and cond2 (cond2 - cond1)
                mean_1:       The mean of cond 1
                mean_2:       The mean of cond 2
                cond_1:       The name of condition1
                cond_2:       The name of condition2
                components:   The components to use (If None, all are used)
                ax:           the axis plotted on
                xrange:       The range on which the axis was plotted on
        """
        
        if components is None:
            components = list(range(len(self.pdata.columns.levels[0])))
        #fi
        components     = np.array(components, dtype='int')
        not_components = np.array([ c for c in range(len(self.pdata.columns.levels[0])) if c not in components ], dtype='int') 

        center_1 = self.pdata.swaplevel(0,1, axis=1).mean()[cond1].values
        center_2 = self.pdata.swaplevel(0,1, axis=1).mean()[cond2].values
        vector = center_2 - center_1
        x = components[0]
        y = components[1]
        slope_2d  = vector[y] / vector[x]
        intercept_2d = ((center_1[y] + center_2[y]) - slope_2d*(center_1[x] + center_2[x])) / 2
        mean_1_2d = ( center_1[x], center_1[y] )
        mean_2_2d = ( center_2[x], center_2[y] )

        vector[not_components]   = 0
        center_1[not_components] = 0
        center_2[not_components] = 0

        xrange = None
        if ax is not None:
            dat    = self.pdata[self.pdata.columns.levels[0][x]][[cond1, cond2]]
            xrange = (dat.min().min(), dat.max().max()) 

            ax.scatter([center_1[x], center_2[x]], [center_1[y], center_2[y]], c=c, zorder=2, label='Center')
            ax.plot(xrange, slope_2d * np.array(xrange) + intercept_2d, c=c, zorder=2, label='Axis of effect')
            ax.legend()
        #fi

        nt = namedtuple('trajectory_axis_of_effect', ['slope_2d', 'intercept_2d', 'mean_1_2d', 'mean_2_2d',
                                                      'vector', 'mean_1', 'mean_2', 'cond_1', 'cond_2',
                                                      'components', 'ax', 'xrange'])

        return nt(slope_2d, intercept_2d, mean_1_2d, mean_2_2d, vector, center_1, center_2,
                  cond1, cond2, components, ax, xrange)
    #edef

#eclass