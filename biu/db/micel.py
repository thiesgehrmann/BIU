from ..structures import Dataset2
from .. import formats
from .. import utils
from .. import ops
from .. import analysis
from ..R import R as biuR

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")
sns = utils.py.loadExternalModule("seaborn")
plt = utils.py.loadExternalModule('matplotlib.pylab')
scipy = utils.py.loadExternalModule('scipy')
sklearn = utils.py.loadExternalModule('sklearn')

import datetime

###############################################################################

class MICEL(Dataset2):
    """
    An interface to the MICEL dataset.
    
    Typical usage:
    --------------
    V = biu.db.MICEL()
    """
    
    def __init__(self, *pargs, **kwargs):

        super(MICEL, self).__init__("micel/", *pargs, **kwargs)
        
        files = [
            ('_d', 'D.csv',  '/home/thies/repos/UA_keelspray/analysis_data/Final_MiCel_Diary.csv'),
            ('_d1', 'D1.csv', '/home/thies/repos/UA_keelspray/analysis_data/MICEL_study_first_survey_index2021_05_05.csv')
        ]
        
        for objname, fname, uri in files:
            self._obj.add_file(fname, utils.Acquire2(uri), finalize=False)
            self._obj.register(objname, [fname], lambda x,d=fname: pd.read_csv(x[d], sep=','))
        #efor
        
        self._obj.add_file('16S.rda', utils.Acquire2('/home/thies/repos/UA_keelspray/analysis_data/MICEL_16S_2021-07-28.rda'), finalize=False)
        self._obj.register('_a', ['16S.rda'], lambda x: biuR()("""load("%s");MICEL$abundance""" % x['16S.rda']))
        self._obj.register('_s', ['16S.rda'], lambda x: biuR()("""load("%s");MICEL$samples""" % x['16S.rda']))
        self._obj.register('T', ['16S.rda'], lambda x: biuR()("""load("%s");MICEL$taxa""" % x['16S.rda']))
        
        def make_s(f):
            S = self._s
            A = self._a
            #S['sample'] = S['sample'].apply(lambda x: x.lower())
            S = S[S.condition.isin(['COV','FAM'])]
            S = S[S['sample_id'].isin(set(A.sample_id) & set(S['sample_id']))]
            S['Participant_number'] = S.apply(lambda x: '%s%03d' % (x.condition, int(x.participant)), axis=1)

            return S.set_index('sample_id')
        #edef
        self._obj.register('S', ['16S.rda'], lambda x: make_s(x['16S.rda']))
        
        def make_a(f):
            S = self.S
            A = self._a
            A = A[A.sample_id.isin(set(A.sample_id) & set(S.index))]
            A = A.pivot(index='sample_id',columns='taxon_id', values='abundance').fillna(0)
            return A
        #edef
        self._obj.register('A', ['16S.rda'], lambda x: make_a(x['16S.rda']))
        
        def make_ra():
            A = self.A
            RA = analysis.microbiome.process.relative(A)
            return RA
        #edef
        
        self._obj.register('RA', [], lambda x: make_ra())

        rename_R = {
             'RecipientFirstNameRecipient_First_Name' : 'Participant_number',
             'RecipientLastNameRecipient_Last_Name' : 'Day_questionaire',
             'Unnamed: 126' : 'Q2_1.Score_Koorts',
             'Q3.4_1Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Hoesten':         'Q3_1.Score_Hoesten',
             'Q3.4_2Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Keelpijn':        'Q3_2.Score_Keelpijn',
             'Q3.4_3Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Neusverkouden':   'Q3_3.Score_Neusverkouden',
             'Q3.4_4Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Kortademig':      'Q3_4.Score_Kortademig',
             'Q3.4_5Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Hoofdpijn':       'Q3_5.Score_Hoofdpijn',
             'Q3.4_6Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Spierpijn':       'Q3_6.Score_Spierpijn',
             'Q3.4_7Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Koude_rillingen': 'Q3_7.Score_Koude_rillingen',
             'Q3.4_8Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Vermoeidheid':  'Q3_8.Score_Vermoeidheid',
             'Q3.4_9Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Verlies_van_reuk/smaak': 'Q3_9.Score_Verlies_van_reuk/smaak',
             'Q3.4_10Scorelijst:_Hoeveel_last_heeft_u_van_onderstaande_klachten?_(1_zeer_mild_en_5_zeer_ernstig_-_Diarree':        'Q3_10.Score_Diarree',
            }
        
        def make_d1(output):

            D = self._d
            D1 = self._d1


            # Rename columns
            D1 = D1.rename(columns=rename_R)
            D1['probiotics'] = (D1[['Q2.23_1Gebruikt_u_probiotica_of_neemt_u_goede_bacteri�n_in_de_onderstaande_vorm:_-_In_zuiverproducten_zoals_Yakult,_Actimel,_Pur_Natur_yoghurt?',
             'Q2.23_2Gebruikt_u_probiotica_of_neemt_u_goede_bacteri�n_in_de_onderstaande_vorm:_-_In_capsules_zoals_Enterol,_Probactiol_of_andere_supplementen_die_bij_de_apotheek__verkrijgbaar_zijn?',]] == 'Ja'
                               ).any(axis=1).astype('category')

            # Fix the dates
            missing_start_dates = dict([
                ("COV003", "25/02/2021"),
                ("COV004", "27/02/2021"),
                ("COV005", "28/02/2021"),
                ("COV008", "04/03/2021"),
                ("COV009", "05/03/2021"),
                ("COV010", "08/03/2021"),
                ("COV013", "11/03/2021"),
                ("COV017", "18/03/2021"),
                ("COV024", "23/03/2021"),
                ("COV015", "16/03/2021"),
                ("COV016", "16/03/2021"),
                ("COV001", "23/02/2021"),
                ("COV006", "01/03/2021"),
                ("COV012", "11/03/2021"),
                ("COV020", "18/03/2021"),
                ("COV021", "22/03/2021")])

            #D1['Q2.15_4_Datum_test'] = D1.apply(lambda x: x['Q2.15_4_Datum_test'] if x['Q2.15_4_Datum_test'] is not None else missing_start_dates[x['Participant_number']],  axis=1)
            D1['Q2.15_4_Datum_test'] = D1.apply(lambda x: missing_start_dates[x['Participant_number']] if x['Participant_number'] in missing_start_dates else x['Q2.15_4_Datum_test'],  axis=1)


            D1['Q2.14_4_Datum_Klachten'] = D1['Q2.14_4_Datum_Klachten'].apply(
                lambda x: datetime.datetime.strptime(x.replace(' ','/'), '%d/%m/%Y') if isinstance(x, str) else x)
            D1['StartDateStart_Date']    = D1['StartDateStart_Date'].apply(
                lambda x: datetime.datetime.strptime(x.split()[0], '%d/%m/%Y'))
            D1['Q2.15_4_Datum_test'] = D1['Q2.15_4_Datum_test'].apply(
                lambda x: datetime.datetime.strptime(x, '%d/%m/%Y') if isinstance(x, str) else x)


            D1['days_to_start'] = (D1['StartDateStart_Date'] - D1['Q2.14_4_Datum_Klachten']).apply(lambda x: x.days)
            D1['days_to_start_from_test'] = (D1['StartDateStart_Date'] - D1['Q2.15_4_Datum_test']).apply(lambda x: x.days)

            D1['days_to_start_missing'] = D1['days_to_start'].isna()
            D1['days_to_start_from_test_missing'] = D1['days_to_start_from_test'].isna()

            D1['BMI'] = D1.apply(lambda x: x['Q2.12Wat_is_uw_gewicht_in_kilogram?']/(x['Q2.11Wat_is_uw_lengte_in_meter___(bijvoorbeeld_1,78)?']**2), axis=1)

            
            R = dict(D.drop_duplicates('Participant_number').set_index('Participant_number')['Randomisation'])
            D1['group'] = [ R.get(i, None) for i in D1.Participant_number ]
            
            d1rename = {
         'Q2.17_1Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Hartaandoening' : 'Hartaandoening',
         'Q2.17_2Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Longaandoening_(bv_astma_of_COPD)': 'Longaandoening',
         'Q2.17_3Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Immuunstoornis': 'Immuunstoornis',
         'Q2.17_4Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Diabetes': 'Diabetes',
         'Q2.17_5Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Reuma': 'Reuma',
         'Q2.17_6Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Hoge_bloeddruk': 'Hoge_bloeddruk',
         'Q2.17_8Heeft_u_een_van_de_volgende_dokter_gediagnostiseerde_chronische_aandoeningen_(meerdere_antwoorden_mogelijk)?_-_Zwaarlijvigheid/obesitas': 'Zwaarlijvigheid/obesitas',
         'Q2.20_1Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Allergische_astma' : 'Allergische_astma',
         'Q2.20_2Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Hooikoorts' : 'Hooikoorts',
         'Q2.20_3Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Huisstofmijt_allergie' : 'Huisstofmijt_allergie',
         'Q2.20_4Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Andere_inhalatieallergie_(huisdieren,_schimmel)' : 'Andere_inhalatieallergie',
         'Q2.20_5Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Eczeem' : 'Eczeem',
         'Q2.20_6Heeft_u_last_van_��n_van_onderstaande_diagnoses_in_relatie_tot_allergie?_Meerdere_antwoorden_zijn_mogelijk_-_Selected_Choice_-_Voedselallergie' : 'Voedselallergie'
            }
            
            D1 = D1.rename(columns=d1rename)
            for c in d1rename.values():
                D1[c] = ~D1[c].isna()
            #efor
            
            
            D1['comorbidity'] = D1[['Hartaandoening', 'Longaandoening', 'Immuunstoornis', 'Diabetes',
                    'Hoge_bloeddruk', 'Zwaarlijvigheid/obesitas', 'Allergische_astma']].any(axis=1)
            
            D1.to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef

        self._obj.add_file('D1.pkl', utils.Acquire2().code(make_d1))
        self._obj.register('D1',['D1.pkl'], lambda x: pd.read_pickle(x['D1.pkl']))

        def make_da(output):

            # Merge all days
            D = self._d
            D1 = self.D1
            DA = pd.concat([D1[rename_R.values()],D[rename_R.values()]])

            DA['Q2_1.Score_Koorts_orig'] = DA['Q2_1.Score_Koorts']
            DA['Q2_1.Score_Koorts'] = DA['Q2_1.Score_Koorts'].apply(lambda x: None if pd.isna(x) else 5 if x > 0 else 0)

            DA['TotalScore'] = DA[['Q2_1.Score_Koorts',
                   'Q3_1.Score_Hoesten', 'Q3_2.Score_Keelpijn', 'Q3_3.Score_Neusverkouden',
                   'Q3_4.Score_Kortademig', 'Q3_5.Score_Hoofdpijn', 'Q3_6.Score_Spierpijn',
                   'Q3_7.Score_Koude_rillingen', 'Q3_8.Score_Vermoeidheid',
                   'Q3_9.Score_Verlies_van_reuk/smaak', 'Q3_10.Score_Diarree']].sum(axis=1)

            DA['SystemScore'] = DA[['Q2_1.Score_Koorts',
                   'Q3_4.Score_Kortademig', 'Q3_6.Score_Spierpijn',
                   'Q3_7.Score_Koude_rillingen', 'Q3_8.Score_Vermoeidheid',
                   'Q3_10.Score_Diarree']].sum(axis=1)

            DA['URTScore'] = DA[[
                   'Q3_1.Score_Hoesten', 'Q3_2.Score_Keelpijn', 'Q3_3.Score_Neusverkouden']].sum(axis=1)

            DA['AcuteScore'] = DA[['Q2_1.Score_Koorts', 'Q3_6.Score_Spierpijn',
                                   'Q3_7.Score_Koude_rillingen', 'Q3_10.Score_Diarree']].sum(axis=1)

            DA = DA[DA.Participant_number.isin(D.Participant_number[~D.Dropout])]
            R = dict(D.drop_duplicates('Participant_number').set_index('Participant_number')['Randomisation'])
            DA['group'] = [ R[i] for i in DA.Participant_number ]

            # Add specific covariate information to the dataframe
            DA = pd.merge(DA, D1[['Participant_number', 
                'Q2.14_4_Datum_Klachten',
                'Q2.15_4_Datum_test',
                'StartDateStart_Date',
                'Q2.1Geslacht_-_Selected_Choice',
                'Q2.3Rookt_u?_-_Selected_Choice',
                'probiotics',
                'Q2.10Wat_is_de_hoogste_opleiding_die_u_met_een_diploma_heeft_afgesloten?',
                'Q2.11Wat_is_uw_lengte_in_meter___(bijvoorbeeld_1,78)?',
                'Q2.12Wat_is_uw_gewicht_in_kilogram?',
                'days_to_start', 'days_to_start_from_test',
                'days_to_start_missing', 'days_to_start_from_test_missing', 'BMI']], on='Participant_number', how='left')

            DA = pd.merge(DA, D1[[ 'Participant_number', 'Hartaandoening', 'Longaandoening', 'Immuunstoornis', 'Diabetes',
                                   'Reuma', 'Hoge_bloeddruk', 'Zwaarlijvigheid/obesitas', 'Allergische_astma', 'Hooikoorts',
                                   'Huisstofmijt_allergie', 'Andere_inhalatieallergie', 'Eczeem', 'Voedselallergie', 'comorbidity']
                                ].set_index('Participant_number').applymap(lambda x: False if pd.isna(x) else True).reset_index(), 
                         on='Participant_number', how='left')


            for c in [ 'Hartaandoening', 'Longaandoening', 'Immuunstoornis', 'Diabetes',
                        'Reuma', 'Hoge_bloeddruk', 'Zwaarlijvigheid/obesitas',
                        'Allergische_astma', 'Hooikoorts', 'Huisstofmijt_allergie',
                        'Andere_inhalatieallergie', 'Eczeem', 'Voedselallergie', 'comorbidity',
                        'days_to_start_missing', 'days_to_start_from_test_missing', ]:
                DA[c] = DA[c].astype('category')
            #efor

            DA['timepoint'] = DA['Day_questionaire'].apply(lambda x: x.replace('D','T'))
            DA.to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef

        self._obj.add_file('DA.pkl', utils.Acquire2().code(make_da))
        self._obj.register('DA',['DA.pkl'], lambda x: pd.read_pickle(x['DA.pkl']))
        
        def make_tsne_embedding(output):
            RA = self.RA
            RA_bc = pd.DataFrame(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(RA, metric='braycurtis')),
                                 index=RA.index, columns=RA.index)

            TRA_bc = pd.DataFrame(sklearn.manifold.TSNE(metric='precomputed').fit_transform(RA_bc),  columns=['D1','D2'], index=RA.index)
            
            TRA_bc.reset_index().set_index('sample_id').to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
        self._obj.add_file('TRA_bc.pkl', utils.Acquire2().code(make_tsne_embedding))
        self._obj.register('TRA_bc',['TRA_bc.pkl'], lambda x: pd.read_pickle(x['TRA_bc.pkl']))
        
        def make_pca_embedding(output):
            RA = self.RA
            RA_bc = pd.DataFrame(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(RA, metric='braycurtis')),
                                 index=RA.index, columns=RA.index)

            PRA_bc = pd.DataFrame(ops.dataframe.pca(RA_bc, nc=2).values, columns=['D1','D2'], index=RA.index)
            
            PRA_bc.reset_index().set_index('sample_id').to_pickle(output)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
        self._obj.add_file('PRA_bc.pkl', utils.Acquire2().code(make_pca_embedding))
        self._obj.register('PRA_bc',['PRA_bc.pkl'], lambda x: pd.read_pickle(x['PRA_bc.pkl']))
    #edef
    
    def plot_embedding(self, hue=None, palette='viridis', data='TSNE',
                       meta=None,  samples=None, s=2,
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
            meta = self.S
        #fi

        if isinstance(data, str):
            if data.upper() == 'TSNE':
                data = self.TRA_bc
            elif data.upper() == 'PCA':
                data = self.PRA_bc
            else:
                data = self.TRA_bc
            #fi
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
        
        if samples is None:
            samples = re.index
        #fi
        re = re.loc[samples]
        
        if hue is not None:
            if re[hue].dtype.name in ['category']:
                re[hue] = re[hue].cat.remove_unused_categories()
            #fi
        #fi
        print(re.shape)

        if ax is None:
            fig, axes = utils.figure.subplots()
            ax = axes[0]
        #fi

        ax = sns.scatterplot(x=d1, y=d2, palette=palette, legend=True, hue=hue, data=re, s=s, ax=ax, *pargs, **kwargs)

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
    
    def selection_strategy(self, cov=True, fam=False):
        DA = self.DA.copy()
        D1 = self.D1.copy()

        print("Total COVID participants:", len(D1.Participant_number.unique()))

        non_dropouts = DA[DA.Day_questionaire.isin(['D21','D20','D19'])].Participant_number.unique()

        DA = DA[DA.Participant_number.isin(non_dropouts)]
        D1 = D1[D1.Participant_number.isin(non_dropouts)]
        #D  = D[D.Participant_number.isin(D.Participant_number[~D.Dropout])]

        print("After removing Dropouts:", len(DA.Participant_number.unique()))

        # Remove people that took more than 6 days to start the study
        DA = DA[(DA.days_to_start_from_test < 4) | (DA.days_to_start < 7)]
        D1 = D1[D1['Participant_number'].isin(DA.Participant_number)]

        # Remove that one bastard with more than 100 days:
        DA = DA[~DA.Participant_number.isin(['COV003'])]
        D1 = D1[~D1.Participant_number.isin(['COV003'])]

        print("After removing those (>=4 days between test and start) and (>=7 days between first symptom and start):", len(DA.Participant_number.unique()))

        # Remove people with a high fever
        DA = DA[~DA.Participant_number.isin(
                 DA[(DA['Day_questionaire'] == 'D1') & (DA['Q2_1.Score_Koorts_orig'] >= 2)].Participant_number.unique())]
        D1 = D1[D1['Participant_number'].isin(DA.Participant_number)]

        print("After removing those with high fever:", len(DA.Participant_number.unique()))

        # Remove people with a too high score at any point (>=40)
        DA = DA[~DA.Participant_number.isin(DA[DA.TotalScore > 40].Participant_number)]
        D1 = D1[D1['Participant_number'].isin(DA.Participant_number)]

        print("After removing those with too high score:", len(DA.Participant_number.unique()))

        participants = list(DA.Participant_number.unique())
        participants = [ p for p in participants if (cov and 'COV' in p) or (fam and 'FAM' in p)]
        return participants
    #edef
#eclass

###############################################################################