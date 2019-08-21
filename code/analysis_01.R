library(here)
library(tidyverse)
library(pROC)

# import data
load(here("data","Full_data.RData"))

# subset clinical into two data frame with patients in the same order as they appear in 
# the gene expression data 
clinicalio <- clinical[clinical$PATIENT_ID %in% io360$PATIENT_ID,]
clinicalio <- clinicalio[match(io360$PATIENT_ID,clinicalio$PATIENT_ID),]

clinicalp <- clinical[clinical$PATIENT_ID %in% pathways$PATIENT_ID,]
clinicalp <- clinicalp[match(pathways$PATIENT_ID,clinicalp$PATIENT_ID),]

#####-----------------------------------------------------------------------------------------#####
#                                     DIFFERENTIAL EXPRESSION
#####-----------------------------------------------------------------------------------------#####

# correlation p-values estimates between levels of histology and each gene ----

# 0 = no neoplastic appearance ; 1 = benign ; 2 = indeterminant ;
# 3 = low-risk malignant ; 4 = high-risk malignant

## CORRELATION ESTIMATES
### IO360 panel 
ioCorr <- apply(io360[,2:ncol(io360)],2,
                function(x)  cor.test(x,clinicalio$histology)$estimate)
### Pathways panel
pathCorr <- apply(pathways[,2:ncol(pathways)],2,
                  function(x)  cor.test(x,clinicalp$histology)$estimate)
## P-VALUES
### IO360 panel
ioCorrP <- apply(io360[,2:ncol(io360)],2,
                 function(x)  cor.test(x,clinicalio$histology)$p.value)
### Pathways panel
pathCorrP <- apply(pathways[,2:ncol(pathways)],2,
                   function(x)  cor.test(x,clinicalp$histology)$p.value)

# combining information into dataset to be used later
## IO360 panel
ioResults <- data.frame(gene = colnames(io360)[2:ncol(io360)],
                  estimate = ioCorr,
                  pval = ioCorrP)
ioResults$gene <- as.character(ioResults$gene)

## Pathways panel 
pathResults <- data.frame(gene = colnames(pathways)[2:ncol(pathways)],
                          estimate = pathCorr,
                          pval = pathCorrP)
pathResults$gene <- as.character(pathResults$gene)

# remove individual pieces calculated
rm(ioCorr, ioCorrP, pathCorr, pathCorrP)

# adjusting for multiple testing ----------------------------------------------
# fdr - least stringent 

## IO360 panel
ioResults$adj_pval <- p.adjust(ioResults$pval,"fdr")
## Pathways panel 
pathResults$adj_pval <- p.adjust(pathResults$pval, "fdr")

# quick summary 
print(paste0("In the IO360 panel ",
             nrow(ioResults[round(ioResults$pval,2)<0.05,]),
             " genes were significantly different (p<0.05) between histologies. Of those, ",
             nrow(ioResults[round(ioResults$adj_pval,2)<0.20,]),
             " remained significant (p<0.2) after adjusting for multiple testing"))
print(paste0("In the Pathways panel ",
             nrow(pathResults[round(pathResults$pval,2)<0.05,]),
             " genes were significantly different (p<0.05) between histologies. Of those, ",
             nrow(pathResults[round(pathResults$adj_pval,2)<0.20,]),
             " remained significant (p<0.2) after adjusting for multiple testing"))

save(ioResults, pathResults,
     file=here("data","Results_data.RData"))