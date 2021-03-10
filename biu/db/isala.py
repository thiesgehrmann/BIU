from ..structures import Dataset2
from .. import formats
from .. import utils
from .. import ops
from .. import R

pd = utils.py.loadExternalModule("pandas")

def _load_all_q(excel):
    X = formats.XLSX(excel, data_only=True)
    Q = X['Sheet1']
    Q.columns = Q.iloc[0]
    Q = Q.iloc[1:]
    Q = Q.set_index(Q.columns[0])

    def try_cast_float(val):
        try:
            return float(val)
        except (ValueError, TypeError) as e:
            return None
        #etry
    #edef

    def try_cast_int(val):
        try:
            return int(val)
        except (ValueError, TypeError) as e:
            return None
        #etry
    #edef

    def try_cast_bool(val):
        if (val is None) or (str(val) == None):
            return None
        else:
            if isinstance(val, bool):
                return val
            elif isinstance(val, str):
                return val.lower() == 'true'
            else:
                return False
            #fi
        #fi
    #edef

    float_cols = [
        'Q1.BMI',
        'Q2.BMI'
    ]
    int_cols = [
        'Q0.Leeftijd',
        'A.day_of_cycle',
        'Q1.Q50.Number_of_pregnancies',
        'Q1.Q6.Number_Of_Biological_Children'
    ]     

    bool_cols = [
        'Q2.Q25.Pill_progestogen',
        'Q2.Q25.Pill_combined',   
        'Q2.Q25.Pill_unknown',
        'Q2.Q25.Spiral_copper',
        'Q2.Q25.Spiral_hormone',
        'Q2.Q25.Prikpil',
        'Q2.Q25.Ring',
        'Q2.Q25.Implant',
        'Q2.Q25.Plaster',
        'Q2.Q25.Condom',
        'Q2.Q25.Cup',
        'Q2.Q25.Morning_after',
        'Q2.Q25.Periodic_celibacy',
        'Q2.Q25.Pullout',
        'Q2.Q25.Sterilized',    
        'Q2.Q25.Partner_sterilized',
        'Q1.Q32.Born_vaginal',
        'Q1.Q32.Born_csect',
        'Q1.Q32.Born_unknown',
        'A.Intercourse_yes',
        'A.Intercourse_no',
        'A.Intercourse_unknown',
    ]

    for col in Q.columns:
        if col in float_cols:
            Q[col] = Q[col].apply(try_cast_float).astype(float)
        elif col in int_cols:
            Q[col] = Q[col].apply(try_cast_int)
        elif col in bool_cols:
            Q[col] = Q[col].apply(try_cast_bool).astype('category')
        #fi
    #efor

    Q['Q2.Delta_BMI']        = Q.apply(lambda x: x['Q2.BMI'] - x['Q1.BMI'], axis=1)
    Q['Q1.Q50.Was_pregnant'] = Q['Q1.Q50.Number_of_pregnancies'].apply(lambda x: True if x > 0 else False)
    Q['Q1.Q6.Has_children']  = Q['Q1.Q6.Number_Of_Biological_Children'].apply(lambda x: True if x > 0 else False)
    Q['A.Critical_period']   = Q['A.day_of_cycle'].apply(lambda x: True if ((x >=14) & (x <= 18)) else False)

    return Q
#edef

class ISALA(Dataset2):
    """
    An interface to the ISALA dataset.
    
    Typical usage:
    --------------
    goto = biu.db.GOTO2()
    
    # RNA-Seq expression data
    goto.expr   # Expression counts
    goto.blood  # Covariates for blood samples
    goto.muscle # Covariates for muscle samples
    
    # Genetic variation
    goto.filter(15, 50001000, 50005000)
    goto.filter_regions([(15, 50001000, 50005000),(13, 50001000, 50005000)])
    
    # Metabolomics data
    goto.metabolomics # Metabolyte measurements
    goto.metabolomics_cov # Covariates for metabolyte data
    """

    data_locations = {
        'Q'     : (_load_all_q,      utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/cleaned_surveys/flow1_all_questionnaires', 'xlsx')),
        'HQRA'  : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/hq_relative_abundances', 'pkl')),
        'HQRAG' : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/hq_relative_abundances_genus', 'pkl')),
        
        'HQS'   : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/hq_samples', 'pkl')),
        
        'PRA'   : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.pca.hqra',  'pkl')),
        'PRAG'  : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.pca.hqrag',  'pkl')),
        
        'TRA'   : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.tsne.hqra',  'pkl')),
        'TRAG'  : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.tsne.hqrag',  'pkl')),
        
        'URA'   : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.umap.hqra',  'pkl')),
        'URAG'  : (pd.read_pickle,   utils.fs.most_recent_file('/home/thies/repos/UA_isala/processed/embedding.umap.hqrag',  'pkl')),
        
        'rds_D'  : (None, '/home/thies/repos/UA_isala/data/16S/isala_cross_amplicon_20210225.rds'),
        'rds_DG' : (None, '/home/thies/repos/UA_isala/data/16S/isala_cross_amplicon_20210225_genus.rds'),
    }
    
    def __init__(self, *pargs, **kwargs):

        super(ISALA, self).__init__("ISALA/", *pargs, **kwargs)
        
        for d in self.data_locations:
            load_func, uri = self.data_locations[d]
            self._obj.add_file(d, utils.Acquire2(uri), finalize=False)
            if load_func is not None:
                self._obj.register(d, [d], lambda x, f=load_func, d=d: f(x[d]))
            #fi
        #efor
        biur = R.R()
        for objname, file, value in [ ('rds_S',  'rds_D', 'samples'),
                             ('rds_A',  'rds_D', 'abundances'),
                             ('rds_T',  'rds_D', 'taxa'),
                             ('rds_AG', 'rds_DG', 'abundances'),
                             ('rds_TG', 'rds_DG', 'taxa') ]:
            self._obj.register(objname, [file], 
                               lambda x, f=file, v=value : biur("""
                                   D <- readRDS('{file:s}')
                                   isala_r_get <- D${value:s}
            """.format(file=x[f], value=v)))
    #edef
#eclass