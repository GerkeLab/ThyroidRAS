# ThyroidRAS
Data and code for examining gene expression differences in RAS thyroid tumors 

## Samples

57 RAS (H,K or N) mutant tumors were surgically resected at Moffitt Cancer Center in Tampa, FL and evaluated using Nanostring's nCounter PanCancer Pathways Panel and IO360 Panel. 

## Nanostring Panels 

The Pathways Panel examines gene expression in 770 genes in 13 cancer related pathways, while the IO360 also examines 770 genes involved in the "interplay between the tumor, microenvironment and immune response".

For more information on these panels, please visit the [Nanostring website](https://www.nanostring.com/products/gene-expression-panels/gene-expression-panels-overview). 

## QC

Raw gene expression was normalized using NanoString nSolver Analysis Software v4.0 and log-2 transformed prior to analysis. 

## Data

The Full_data.RData file contains 5 text files also availble in the data folder. Clinical.txt contains clinical features associated with each tumor, while Pathways.txt and IO360.txt contains log2 expression values for all genes on each panel (some genes overlap between panels). The PathwaysAnnotation.txt and IO360Annotation.txt files denote which genes are in which pathways using +/- annotation. 

## Code 

The analysis is broken into two scripts analysis_01 which calculates differential expression and the scaled risk score and figures_02 which uses the results from analysis_01 to reproduce the figures from the manuscript. Results from analysis_01 can also be founf in the Results_data.RData file. 

## Collaborators

Juan C. Hernandez-Prera
Pablo Valderrabano
Jordan Creed
Barbara Centeno
Valentina Tarasova
Julie Hallanger-Johnson
Bryan McIver
Bruce Wenig
Sean Yoder
Cesar A. Lam
Derek S. Park
Alexander R. Anderson
Natarajan Raghunand
Anders Berglund
Travis A. Gerke
Christine H. Chung


## Contact Information 

Questions or comments about the data or code in this repo can be directed to Travis Gerke @ travis.gerke@moffitt.org or Jordan Creed @ Jordan.H.Creed@moffitt.org .
