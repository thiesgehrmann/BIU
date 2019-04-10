from ... import utils
from ... import ops

sklearn = utils.py.loadExternalModule('sklearn')

#################################################################################

def impute_find_train(df, col):
    """
    Finds the training set to use for a specific column.
    Rule: for each row with nan value in the column, we must find other columns that are not nan in the rows for which there are no nans in the specified column.
    """
    col_filled = df[~df[col].isna()]
    col_nan    = df[df[col].isna()]
    
    # Identify columns that are not nan in ALL col_nan rows,
    # and not nan in ALL col_filled rows.
    
    nnc_filled = col_filled.columns[~col_filled.isna().any(0)]
    nnc_nan    = col_nan.columns[~col_nan.isna().any(0)]
    
    train = col_filled[list(set(list(nnc_filled & nnc_nan) + [col]))]

    return train
#edef

#################################################################################

def impute_numeric_nn(train, pred, col):
    """
    A nearest neighbour imputation scheme.
    """
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(train.drop(columns=[col]))
    
    nnbrs = ops.lst.flatten(nbrs.kneighbors(pred.drop(columns=[col]).values, 1, return_distance=False))
    return train[col].values[nnbrs]
#edef

#################################################################################
    

def impute_numeric(df, columns=None, maxiter=10, imputer='nn'):
    """
    Perform an imputation on numerical data
    
    parameters:
    -----------
    df: Pandas DataFrame. The dataframe to perform the imputation on
    columns: The columns to impute (if None, then all columns with Nan values will be imputed)
    maxiter: The maximum number of iterations to perform in the imputation.
    imputer: String|Function. The imputer to use. Must be in ['nn'], or an object with __call__ attribute.
    
    returns:
    Pandas Dataframe with imputed columns
    """
    
    defined_imputers = {
            'nn' : impute_numeric_nn
        }
    
    if not hasattr(imputer, '__call__'):
        imputer = imputer.lower()
        if imputer not in defined_imputers:
            raise ValueError("'%s' is not a defined imputer." % imputer)
        #fi
        imputer = defined_imputers[imputer]
    #fi
    
    if columns is None:
        columns = df.columns[df.isna().any(0)]
    #fi
    
    imputed = df.copy()
    last_imputed = df.copy()
    
    for i in range(maxiter):
        
        for col in columns:
            train_data = impute_find_train(imputed, col)
            pred_data  = imputed.loc[df[col].isna()][train_data.columns]

            imputed[col][df[col].isna()] = imputer(train_data, pred_data, col)

        #efor
        
        if imputed.equals(last_imputed):
            break
        #fi
        
        last_imputed = imputed.copy()
    #efor
    
    return imputed
#edef

#################################################################################

def impute_categorical(df, columns=None, maxiter=10, impute_func=None):
    raise NotImplementedError("This functionality is not yet implemented...")
#edef

#################################################################################
