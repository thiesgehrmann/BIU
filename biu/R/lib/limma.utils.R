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
    
    res <- if (!('contrasts' %in% names(fit))){
            r <- data.frame(topTreat(fit, coef=contr.idx, n=Inf))
            r$contr <- "main effect"
            r$gene <- rownames(r)
            r
        } else {
    
            if (is.null(contr.names)){
                contr.names <- colnames(fit$contrasts)
            }

            all_test_res = NULL
            for (contr.idx in 1:dim(fit$contrasts)[2]){
                contr.name <- contr.names[contr.idx]
                testRes <- data.frame(topTreat(fit, coef=contr.idx, n=Inf))
                testRes$contr <- contr.name
                testRes$gene  <- rownames(testRes)
                if (is.null(all_test_res)){
                    all_test_res <- testRes
                } else {
                    all_test_res <- rbind(all_test_res, testRes)
                }
            }
            all_test_res
    
        }

    colnames(res)[colnames(res)=="P.Value"] <- "pvalue"
    colnames(res)[colnames(res)=="adj.P.Val"] <- "qvalue"
    res$fdr <- p.adjust(res$pvalue, method='fdr')
    
    res
    }

#########################################################################

limma.metab <- function(formula, covar, metab, contrasts=NULL, random_effect=NULL){
    #'
    #'parameters:
    #'-----------
    #'formula: The formula to evaluate
    #'covar: DataFrame Covariates (defined in covar)
    #'metab: Dataframe of Metabolite measurements (Normalized)
    #'contrasts: list. Which contrasts to evaluate (default NULL)
    #'random_effect: String. Add a random effect per group defined by column in covar 
    
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

    limma.extract.contrasts(fit)
}

#########################################################################

limma.rnaseq <- function(formula, covar, expr, contrasts=NULL, random_effect=NULL, voom=TRUE){
    #'
    #'parameters:
    #'-----------
    #'formula: The formula to evaluate
    #'covar: DataFrame Covariates (defined in covar)
    #'expr: Dataframe of count measurements
    #'contrasts: list. Which contrasts to evaluate (default NULL)
    #'random_effect: String. Add a random effect per group defined by column in covar 
    
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

    # Add covariates
    for (cov in names(C)){
      dge$samples[[cov]] <- C[[cov]]
    }
    
    design <- model.matrix(formula, dge$samples)
    
    dge <- calcNormFactors(dge)
    if (voom){
        dge <- voom(dge, design)
    }
    
    dup.corr <- if (!is.null(random_effect)){
                    print("Estimating correlation")
                    d <- duplicateCorrelation(dge, design=design, block=dge$samples[,random_effect])
                    print(d$consensus)
                    d
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

    limma.extract.contrasts(fit)
    }