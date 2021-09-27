from ..structures import Dataset2
from .. import formats
from .. import utils
from .. import ops
from .. import analysis

pd = utils.py.loadExternalModule("pandas")
sns = utils.py.loadExternalModule("seaborn")
plt = utils.py.loadExternalModule('matplotlib.pylab')
scipy = utils.py.loadExternalModule('scipy')
sklearn = utils.py.loadExternalModule('sklearn')


###############################################################################

class Valencia(Dataset2):
    """
    An interface to the Ravel VALENCIA dataset.
    
    Typical usage:
    --------------
    V = biu.db.Valencia()
    """

    data_locations = [
        ('_counts' , 'counts.csv', 'https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/Publication_materials/Data_and_metadata/all_samples_taxonomic_composition_data.csv'),
        ('_meta'   , 'meta.csv', 'https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/Publication_materials/Data_and_metadata/all_samples_metadata.csv'),
    ]
        
    
    def __init__(self, *pargs, **kwargs):

        super(Valencia, self).__init__("valencia/", *pargs, **kwargs)
        
        for objname, fname, uri in self.data_locations:
            self._obj.add_file(fname, utils.Acquire2().curl(uri))
            self._obj.register(objname, [fname], lambda x,d=fname :pd.read_csv(x[d], sep=','))
        #efor
        
        def make_RA(output):
            D = self._counts.set_index('Sample_number_for_SRA').drop(
                    columns=['Subject_number','HC_CST','HC_subCST','Val_CST','Val_subCST','total_reads',
                             'I-A_sim', 'I-B_sim', 'II_sim', 'III-A_sim', 'III-B_sim', 'IV-A_sim', 'IV-B_sim',
                             'IV-C0_sim', 'IV-C1_sim', 'IV-C2_sim', 'IV-C3_sim', 'IV-C4_sim', 'V_sim'])
            RA = analysis.microbiome.process.relative(D)
            RA = RA.reset_index().rename(columns={'Sample_number_for_SRA':'index'}).set_index('index')
            RA.to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
        self._obj.add_file('RA.pkl', utils.Acquire2().code(make_RA))
        self._obj.register('RA',['RA.pkl'], lambda x: pd.read_pickle(x['RA.pkl']))

        def make_meta(output):
            M = self._meta.set_index('Sample_number')
            C = self._counts.set_index('Sample_number_for_SRA')[
                ['HC_CST','HC_subCST','Val_CST','Val_subCST',
                 'I-A_sim', 'I-B_sim', 'II_sim', 'III-A_sim', 'III-B_sim', 'IV-A_sim', 'IV-B_sim',
                 'IV-C0_sim', 'IV-C1_sim', 'IV-C2_sim', 'IV-C3_sim', 'IV-C4_sim', 'V_sim']]
            RA = self.RA
            
            J = C.join(M, how='outer').reset_index().drop_duplicates('index').set_index('index')
            J.loc[RA.index, 'ra_genus_1'] = RA.apply(lambda r: RA.columns[ops.lst.argrank(r,0)], axis=1)
            J.loc[RA.index, 'ra_genus_2'] = RA.apply(lambda r: RA.columns[ops.lst.argrank(r,1)], axis=1)
            J.loc[RA.index, 'ra_genus_3'] = RA.apply(lambda r: RA.columns[ops.lst.argrank(r,2)], axis=1)
            J.to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
        self._obj.add_file('RAS.pkl', utils.Acquire2().code(make_meta))
        self._obj.register('RAS',['RAS.pkl'], lambda x: pd.read_pickle(x['RAS.pkl']))
        
        def make_embedding(output):
            RA = self.RA
            RA_bc = pd.DataFrame(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(RA, metric='braycurtis')),
                                 index=RA.index, columns=RA.index)

            TRA_bc = pd.DataFrame(sklearn.manifold.TSNE(metric='precomputed').fit_transform(RA_bc),  columns=['D1','D2'], index=RA.index)
            
            TRA_bc.reset_index().set_index('index').to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
        self._obj.add_file('TRA_bc.pkl', utils.Acquire2().code(make_embedding))
        self._obj.register('TRA_bc',['TRA_bc.pkl'], lambda x: pd.read_pickle(x['TRA_bc.pkl']))

    #edef

    
    def plot_embedding(self, hue=None, palette='viridis', data='TSNE',
                       meta=None,
                       legend=False, title=None, d1='D1', d2='D2', ax=None, *pargs, **kwargs):
        """
        plot_embedding

        parameters:
        -----------
        hue : string|list|None
            Describes colors for each point
            if string, then column in sample_meta or meta
            if list, then list of strings or numbers, of length equal to data
            if None, then no color is plotted
        palette : string
            color palette to use for hue, valid for plt.get_cmap
        data: pd.DataFrame|'TSNE'|'UMAP'|'PCA'
            dataset with embedding points to plot
            if pd.Dataframe, then use this data.
            if 'TSNE', use tsne embedding
            if 'UMAP', use umap embedding
            otherwise, use PCA embedding
        sample_meta : pd.DataFrame
            A dataframe with contains information PER SAMPLE IN DATA
        meta: pd.DataFrame
            A dataframe which contains information per sample, but is indexed on a different column
        meta_match : String
            The column by which meta is indexed, but exists in sample_meta
        legend: bool
            Plot a legend or not
        title: String
            Title for the plot
        d1,d2: Strings
            The names of the first and second components of the embedding to use
        *pargs, **kwargs : dicts
            additional arguments for seaborn.scatterplot
            INCLUDES ax!

        """

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes


        if meta is None:
            meta = self.RAS
        #fi

        if isinstance(data, str):
            data = self.TRA_bc
        #fi

        re = data.copy()
        if hue is not None:
            try:
                if hue in meta.columns:
                    re[hue] = meta[hue].values
                #fi
            except TypeError:
                    re['hue'] = hue
                    hue = 'hue'
                #etry
            #fi
        #fi

        if ax is None:
            fig, axes = utils.figure.subplots()
            ax = axes[0]
        #fi

        ax = sns.scatterplot(x=d1, y=d2, palette=palette, legend=True, hue=hue, data=re, s=2, ax=ax, *pargs, **kwargs)

        if legend:
            if pd.api.types.is_bool_dtype(re[hue]):
                ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol= 2)
            elif pd.api.types.is_numeric_dtype(re[hue]):
                norm = plt.Normalize(re[hue].min(), re[hue].max())
                sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
                if ax.get_legend() is not None:
                    ax.get_legend().remove()
                #fi

                axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
                ax.figure.colorbar(sm, cax=axins)
            else:
                ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol= 2)
            #fi
        elif hue is not None:
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            #fi
        #fi

        if title:
            ax.set_title(title)
        #fi

        return ax
    #edef
#eclass

###############################################################################
