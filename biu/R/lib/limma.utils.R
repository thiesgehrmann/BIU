library(limma)
library(edgeR)

#########################################################################

normalize.rin <- function(x){
    ## Defining a function for Rank Inverse Normal transformation of a single measurement:
    single.RIN <- function(x){
        x <- rank(x, "keep")
        x  <- (x - 0.5) / length(x)
        return(qnorm(x))
    }
    if(is.matrix(x)){
        x <- t(apply(x,1,single.RIN))
        return(x)
    } else if(is.data.frame(x)){
        y <- t(apply(x,1,single.RIN))
        y <- data.frame(y)
        rownames(y) <- rownames(x)
        colnames(y) <- colnames(y)
        return(y)
    } else {
        if(is.vector(x)){
            return(single.RIN(x))
        }
    }
}

#########################################################################

limma.extract.contrasts <- function(fit, contr.names=NULL){
    #' For the output from eBayes, extract all the tests, per contrast that were preent in the object.
    #' parameters:
    #' -----------
    #' fit2C: The output from eBayes
    #' contr.names: A list of alternative names for the contrasts.
    #'
    #' Returns:
    #' --------
    #' A dataframe with test results

    # Combine outputs
    
    res <- NULL
    if (is.null(contr.names)){
        contr.names <- colnames(fit$contrasts)
    }

    for (contr.idx in 1:dim(fit$contrasts)[2]){
        contr.name <- contr.names[contr.idx]
        contr.res <- data.frame(topTreat(fit, coef=contr.idx, n=Inf))
        contr.res$contr <- contr.name
        contr.res$gene  <- rownames(contr.res)
        if (is.null(res)){
            res <- contr.res
        } else {
            res <- rbind(res, contr.res)
        }
    }

    colnames(res)[colnames(res)=="P.Value"] <- "pvalue"
    colnames(res)[colnames(res)=="adj.P.Val"] <- "qvalue"
    res$fdr <- p.adjust(res$pvalue, method='fdr')
    
    res
    }

limma.extract.effects <- function(fit, effects=NULL) {
    #' For the output from eBayes, extract all the tests, per effect defined in the formula (or provided in the effects list)
    #' parameters:
    #' -----------
    #' fit2C: The output from eBayes
    #' effects: A list of effects in which we are interested.
    #'
    #' Returns:
    #' --------
    #' A dataframe with test results

    if (is.null(effects)){
        effects <- colnames(fit$design)
    }
    
    res = NULL
    for (effect.name in effects){
        effect.res <- data.frame(topTreat(fit, coef=effect.name, n=Inf))
        effect.res$contr <- effect.name
        effect.res$gene  <- rownames(effect.res)
        if (is.null(res)){
            res <- effect.res
        } else {
            res <- rbind(res, effect.res)
        }
    }
    
    colnames(res)[colnames(res)=="P.Value"] <- "pvalue"
    colnames(res)[colnames(res)=="adj.P.Val"] <- "qvalue"
    res$fdr <- p.adjust(res$pvalue, method='fdr')
    
    res
    
    }

#########################################################################

limma.metab <- function(formula, covar, metab, contrasts=NULL, effects=NULL, random_effect=NULL, verbose=FALSE){
    #'
    #'parameters:
    #'-----------
    #'formula: The formula to evaluate
    #'covar: DataFrame Covariates (defined in covar)
    #'metab: Dataframe of Metabolite measurements (Normalized)
    #'contrasts: list. Which contrasts to evaluate (default NULL)
    #'effects: list. Which effects are you interested in? (You cannot provide contrasts AND effects. You must choose one).
    #'random_effect: String. Add a random effect per group defined by column in covar 
    #'verbose: Boolean. Return all the intermediate variables in a list

    
    C <- na.omit(droplevels(covar))
    E <- as.matrix(metab[match(rownames(C), rownames(metab)),])
    E <- apply(E,2,scale)
    E <- t(normalize.rin(E))
    design <- model.matrix(formula, C)
    
    dup.corr <- if (!is.null(random_effect)){
                    dup.corr <- duplicateCorrelation(E, design=design, block=C[,random_effect])
                } else{
                    NULL
                }

    contr.matrix <- if (! is.null(contrasts)){
                        cm  <- makeContrasts(contrasts=as.vector(contrasts), levels=design)
                        if (!is.null(names(contrasts))){
                            colnames(cm) <- names(contrasts)
                        }
                        cm
                    } else {
                        NULL
                    }

    # Fit limma
    fit <- if (!is.null(random_effect)){
                lmFit(E, design, block=C[,random_effect], correlation=dup.corr$consensus)
            } else{
                lmFit(E, design)
            }
    fit <- if (!is.null(contr.matrix)){
               fit2C <- contrasts.fit(fit, contr.matrix)
           } else {fit}
    
    fit <- eBayes(fit)

    res <- if (! is.null(contrasts)){
        limma.extract.contrasts(fit)
        } else {
        limma.extract.effects(fit, effects)
        }

    if (verbose){
        return(list(res=res,
                    fit=fit,
                    contr.matrix=contr.matrix,
                    dup.corr=dup.corr))
    } else {
        return(res)
    }
}

#########################################################################

limma.rnaseq <- function(formula, covar, expr, contrasts=NULL, effects=NULL, random_effect=NULL, voom=TRUE, verbose=FALSE){
    #'
    #'parameters:
    #'-----------
    #'formula: The formula to evaluate
    #'covar: DataFrame Covariates (defined in covar)
    #'expr: Dataframe of count measurements
    #'contrasts: list. Which contrasts to evaluate (default NULL)
    #'effects: list. Which effects are you interested in? (You cannot provide contrasts AND effects. You must choose one).
    #'random_effect: String. Add a random effect per group defined by column in covar 
    #'voom: Boolean. Perform VOOM normalization
    #'verbose: Boolean. Return all the intermediate variables in a list
    
    C <- na.omit(droplevels(covar))
    E <- expr[ rownames(expr) %in% rownames(C),]

    # Make sure that the data frames are in the same order!
    C <- C[order(rownames(C)),]
    E <- E[order(rownames(E)),]

    ##############################################################
    # Prepare the Diff. Ex. data structures

    # Put gene expression matrix into DGEList object. Remove the sample column
    dge <- DGEList(t(E))
    # Add the sample names to the DGEList
    colnames(dge) <- rownames(E)

    design <- model.matrix(formula, C)
    
    dge <- calcNormFactors(dge)
    if (voom){
        dge <- voom(dge, design)
    }
    
    dup.corr <- if (!is.null(random_effect)){
                    duplicateCorrelation(dge, design=design, block=C[,random_effect])
                } else{
                    NULL
                }

    contr.matrix <- if (! is.null(contrasts)){
                        cm  <- makeContrasts(contrasts=as.vector(contrasts), levels=design)
                        if (!is.null(names(contrasts))){
                            colnames(cm) <- names(contrasts)
                        }
                        cm
                    } else {
                        NULL
                    }

    # Fit limma
    fit <- if (!is.null(random_effect)){
                lmFit(dge, design, block=C[,random_effect], correlation=dup.corr$consensus)
            } else{
                lmFit(dge, design)
            }
    fit <- if (!is.null(contr.matrix)){
                contrasts.fit(fit, contr.matrix)
           } else {fit}
    
    fit <- eBayes(fit)
    
    res <- if (! is.null(contrasts)){
        limma.extract.contrasts(fit)
        } else {
        limma.extract.effects(fit, effects)
        }

    if (verbose){
        return(list(res=res,
                    fit=fit,
                    contr.matrix=contr.matrix,
                    dup.corr=dup.corr))
    } else {
        return(res)
    }
    
    }