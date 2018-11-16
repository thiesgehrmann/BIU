from ... import utils

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

def framingham_risk_score(male, age, cholesterol, HDL, SBP, HT_treatment, smoking_status, diabetes_status):
    """
    Calculate the cardiovascular FRS for an individual
    Inputs:
        male:            Boolean. The sex of the individual. True if male.
        age:             Float. The age of the invididual
        cholesterol:     Float. Cholesterol levels
        HDL:             Float. High Density LipoProtein levels
        SBP:             Float. Systolic Blood Pressure
        HT_treatment:    Boolean. Does this person take anti-hypertensive medication? True if yes.
        smoking_status:  Boolean. Does this person smoke? True if yes.
        diabetes_status: Boolean. Does this person have diabetes? True if yes.
    Output:
        FRS score.
        Based on Cox model method from https://doi.org/10.1161/CIRCULATIONAHA.107.699579

        returns None if any value is missing or NA
        
    # Case 1 from paper: Expect 0.1048
    framingham_risk_score(False, 61, 180, 47, 124, False, True, False)
    # > 0.10484443957186262

    # Case 2 from paper: Expect 0.1562
    framingham_risk_score(True, 53, 161, 55, 125, True, False, True)
    # > 0.156231838948365
    """

    if np.any(pd.isna([male, age, cholesterol, HDL, SBP, HT_treatment, smoking_status, diabetes_status])):
      return None
    #fi

    male_params    = np.array([ 3.06117, 1.12370, -0.93263, 1.93303, 1.99881, 0.65451, 0.57367 ])
    male_average   = np.array([ 3.8560,  5.3420,   3.7686,  4.3544,  0.5019,  0.3522,  0.0650  ])
    female_params  = np.array([ 2.32888, 1.20904, -0.70833, 2.76157, 2.82263, 0.52873, 0.69154 ])
    female_average = np.array([ 3.8686,  5.3504,   4.0176,  4.2400,  0.5826,  0.3423,  0.0376  ])
    
    HT_treatment    = 1 if HT_treatment    else 0
    smoking_status  = 1 if smoking_status  else 0
    diabetes_status = 1 if diabetes_status else 0

    values = np.array([ np.log(age),
                        np.log(cholesterol),
                        np.log(HDL),
                        np.log(SBP)*(1-HT_treatment),
                        np.log(SBP)*HT_treatment,
                        smoking_status,
                        diabetes_status ])

    if male:
        effect = np.sum((values - male_average) * male_params)
        return 1 - np.power(0.88936, np.exp(effect))
    else:
        effect = np.sum((values - female_average) * female_params)
        return 1 - np.power(0.95012, np.exp(effect))
    #fi
#edef
