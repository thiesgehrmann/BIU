library('tidyr')
#library('lmer')

lm.remove.effects <- function(F, D, L){
    #' Returns the residuals of a model.
    #' parameters:
    #' -----------
    #' F: Formula
    #' D: The data from which you want to remove the effects described in F
    #'    Rows are samples, columns are measurements (e.g. mrna transcripts or metabolites)
    #' L: Covariates, corresponding to the variables in L
    #'
    #' A combined model is learned across all measurements (a covariate, e.g. sex, has the same effect in each metabolite)
    #'

    
    D <- droplevels(D)
    L <- droplevels(L)

    D$lm.remove.effects.row.id <- rownames(D)
    L$lm.remove.effects.row.id <- rownames(L)

    DL  <- merge(D, L, on='lm.remove.effects.row.id')
    PDL <- gather(DL, colnames(D)[! colnames(D) %in% colnames(L)],
                  key='lm.remove.effects.meas',
                  value='lm.remove.effects.level')


    F <- update(F, lm.remove.effects.level ~ . + lm.remove.effects.meas)

    m <- lm(F, PDL)

    PDL$lm.remove.effects.resid <- residuals(m)
    PDL$lm.remove.effects.level <- NULL
    
    DL.spread <- spread(PDL, lm.remove.effects.meas, lm.remove.effects.resid)
    rownames(DL.spread) <- DL.spread$lm.remove.effects.row.id
    DL.spread <-DL.spread[colnames(D)[! colnames(D) %in% colnames(L)]]
    
    DL.spread
    }