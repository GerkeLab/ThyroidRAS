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

# correlation p-values estimates between levels of histology and each gene 
# 0 = no neoplastic appearance ; 1 = benign ; 2 = indeterminant ; 3 = low-risk malignant ;
# 4 = high-risk malignant
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

# adjusting for multiple testing 
# fdr since least stringent 
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

#####-----------------------------------------------------------------------------------------#####
#                                             PATHWAYS
#####-----------------------------------------------------------------------------------------#####

# get pathways names 
io_paths <- colnames(io360_panel)[2:ncol(io360_panel)]
path_paths <- colnames(pathways_panel)[2:ncol(pathways_panel)]

# calculate risk scores by adding expression values for all genes significantly assocaited with 
# histology. Negative weights applied to those inversely associated with increasing agressiveness
## IO360 panel
io_score <- c()
for (i in io_paths){
  gene_names <- ioResults[round(ioResults$pval,2)<0.05 &
                      (ioResults$gene %in% io360_panel[io360_panel[,i]=="+",]$Gene), ]$gene
  pos_assoc <- ioResults[ioResults$estimate > 0,]$gene
  neg_assoc <- ioResults[ioResults$estimate < 0,]$gene
  if(length(gene_names)>1){
    if(length(intersect(gene_names,pos_assoc)) > 1){
      tmp1 <- apply(io360[,intersect(gene_names,pos_assoc)], 1, function(x) sum(x))
    } else if (length(intersect(gene_names,pos_assoc))==1) {
      tmp1 <- io360[,intersect(gene_names,pos_assoc)]
    } else {tmp1 <- 0}
    
    if(length(intersect(gene_names,neg_assoc)) > 1){
      tmp2 <- apply(io360[,intersect(gene_names,neg_assoc)], 1, function(x) sum(x))
    } else if (length(intersect(gene_names,neg_assoc))==1) {
      tmp2 <- io360[,intersect(gene_names,neg_assoc)]
    } else {tmp2 <- 0}
    
    tmp <- tmp1 - tmp2
  } else { tmp <- io360[,gene_names]}
  io_score <- as.data.frame(cbind(io_score,tmp))
}
## pathways panel 
path_score <- c()
for (i in path_paths){
  gene_names <- pathResults[round(pathResults$pval,2)<0.05 &
                        (pathResults$gene %in% pathways_panel[pathways_panel[,i]=="+",]$Gene), ]$gene
  pos_assoc <- pathResults[pathResults$estimate > 0,]$gene
  neg_assoc <- pathResults[pathResults$estimate < 0,]$gene
  if(length(gene_names)>1){
    if(length(intersect(gene_names,pos_assoc)) > 1){
      tmp1 <- apply(pathways[,intersect(gene_names,pos_assoc)], 1, function(x) sum(x))
    } else if (length(intersect(gene_names,pos_assoc))==1) {
      tmp1 <- pathways[,intersect(gene_names,pos_assoc)]
    } else {tmp1 <- 0}
    
    if(length(intersect(gene_names,neg_assoc)) > 1){
      tmp2 <- apply(pathways[,intersect(gene_names,neg_assoc)], 1, function(x) sum(x))
    } else if (length(intersect(gene_names,neg_assoc))==1) {
      tmp2 <- pathways[,intersect(gene_names,neg_assoc)]
    } else {tmp2 <- 0}
    
    tmp <- tmp1 - tmp2
  } else { tmp <- pathways[,gene_names]}
  path_score <- as.data.frame(cbind(path_score,tmp))
}

# remove individual parts from environment
rm(gene_names, pos_assoc, neg_assoc, tmp, tmp1, tmp2, i)

colnames(path_score) <- path_paths # add name for pathways
path_score <- scale(path_score) # scale risk scores - easier to compare
path_score <- data.frame(PATIENT_ID = clinicalp$PATIENT_ID, # data frame to use later 
                         histology = clinicalp$histology,
                         path_score)

colnames(io_score) <- io_paths # add name for pathways
io_score <- scale(io_score) # scale risk scores - easier to compare
io_score <- data.frame(PATIENT_ID = clinicalio$PATIENT_ID, # data frame to use later 
                       histology = clinicalio$histology,
                       io_score)

# examining how well risk scores correlate with histology aggressiveness
# calcualting correlation estimates and pvalues 
# IO360 panel
io_risk_p <- c()
io_risk_corr <- c()
for(i in gsub(" ","\\.",gsub("-","\\.",io_paths))){
  io_risk_corr <-  c(io_risk_corr, cor.test(io_score[,i],io_score$histology)$estimate)
  io_risk_p <- c(io_risk_p, cor.test(io_score[,i],io_score$histology)$p.value)
}
## in pathways panel
path_risk_p <- c()
path_risk_corr <- c()
for(i in gsub(" ","\\.",gsub("-","\\.",gsub("\\+","\\.",path_paths)))){
  path_risk_corr <-  c(path_risk_corr, cor.test(path_score[,i],path_score$histology)$estimate)
  path_risk_p <- c(path_risk_p, cor.test(path_score[,i],path_score$histology)$p.value)
}

io_risk <- data.frame(pathway = gsub(" ","\\.",gsub("-","\\.",io_paths)),
                      estimate = io_risk_corr,
                      pval = io_risk_p)

path_risk <- data.frame(pathway = gsub(" ","\\.",gsub("-","\\.",gsub("\\+","\\.",path_paths))),
                        estimate = path_risk_corr,
                        pval = path_risk_p)

rm(i, path_risk_corr, path_risk_corr, io_risk_corr, io_risk_p)

save(ioResults, pathResults, io_score, path_score, io_risk, path_risk,
     file="/Volumes/Lab_Gerke/thyroidInnovation/github/data/Results_data.RData")