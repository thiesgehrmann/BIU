#!/usr/bin/env Rscript

##############################################################

library('edgeR')
library('limma')

##############################################################


args = commandArgs(trailingOnly=TRUE)

formula   = args[1]
group     = args[2]
contrasts = args[3]
Efile     = args[4]
Cfile     = args[5]
Ofile     = args[6]

##############################################################
# Process inputs

# Read the expression data
E <- read.table(Efile, sep=';', row.names=1, header=TRUE)

# Reads the covariates
C <- read.table(Cfile, sep=';', row.names=1, header=TRUE, stringsAsFactors=TRUE)

# Prepare the formula
F <- as.formula(formula)

# Prepare the contrasts
contrasts <- unlist(strsplit(contrasts, ";"))

##############################################################
# Select the data, etc

# Only select samples for which we have the covariate information
E <- E[ rownames(E) %in% rownames(C),]

# Make sure that the data frames are in the same order!
C <- C[order(rownames(C)),]
E <- E[order(rownames(E)),]

##############################################################
# Prepare the Diff. Ex. data structures

# Put gene expression matrix into DGEList object. Remove the sample column
dge <- DGEList(t(E))

# Add the sample names to the DGEList
colnames(dge) <- rownames(E)

# Add covariates. Skip sample and individual These are not necessary as covariates.
for (cov in names(C)){
  dge$samples[[cov]] <- C[[cov]]
}

##############################################################
# Create the model

# Create the design matrix. Factors are expanded. For our specific 
design <- model.matrix(F, data=dge$samples)
colnames(design) <- gsub(group, "", colnames(design))

# Define the differential expression tests were are interested in
contrast.list = list()
for (item in 1:length(contrasts)){
  contrast.list[[item]] = c(contrasts[item])
}
contrast.list[["levels"]] = colnames(design)

contr.matrix <- do.call(makeContrasts, contrast.list)

##############################################################
# Normalize the expression

dgeTMM  <- calcNormFactors(dge)
dgeVoom <- voom(dgeTMM, design)

##############################################################
# Fit the model

voomfit <- lmFit(dgeVoom, design)
voomfit <- contrasts.fit(voomfit, contrasts=contr.matrix)
efit <- eBayes(voomfit)

##############################################################
# Perform the tests

#plotSA(efit)
#summary(decideTests(efit))

#tfit <- treat(efit, lfc=.5)
#summary(decideTests(tfit))

##############################################################
# Output the data

allTestRes = NULL
for (contr.idx in 1:dim(contr.matrix)[2] ){
  contr.name <- names(contr.matrix[1,])[contr.idx]
  testRes <- data.frame(topTreat(efit, coef=contr.idx, n=Inf))
  testRes$contr <- unlist(strsplit(contr.name,'='))[[1]]
  testRes$gene  <- rownames(testRes)
  print(contr.name)
  
  if (is.null(allTestRes)){
    allTestRes <- testRes
  } else {
    allTestRes <- rbind(allTestRes, testRes)
  }
}

write.table(allTestRes, file=Ofile, sep=';', row.names=FALSE)

##############################################################

