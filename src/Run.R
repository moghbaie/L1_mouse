## Mehrnoosh Oghbaie
## 07/24/2019

## Define class and run the project from this file
rm(list=ls())

setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

CRAN.packages <- c("readr","readxl","data.table","reshape2","dplyr","magrittr",
                   "igraph","sqldf","stringr","corrplot","ggplot2","R6","ggridges",
                   "gridExtra","ggrepel","rgl","venn", 'writexl')
bioconductor.packages <- c("biomaRt","limma","qvalue","msa","ape","seqinr","ggseqlogo","PTXQC")


############################################################################################
## Check the values in input.info file

source("functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("functions/Data_preparation.R")
source("functions/Imputation.R")
source("functions/Anova.R")

run_order <- rbind(cbind("190607_E_FreyaLC_JP_V5","190607_E_FreyaLC_JP_IgG"),
                   cbind("tk180324_JLC_L._V5","tk180324_JLC_L._IgG"), 
                   cbind("tk180324_JLC_S._V5","tk180324_JLC_S._IgG"),
                   cbind("tk180325_JLC_Brain._V5","tk180325_JLC_Brain._IgG"),
                   cbind("tk180325_JLC_Neonate._V5","tk180325_JLC_Neonate._IgG"))

## Data_preparation.R
L1 <- Template$new()
## Produce quality control report in input directory - reads it from input.info
#runQC(L1$input)

L1$removeContaminant()
L1$logTransformation()
L1$separatedList(run_order)

## Imputation.R
L1$removeAllZeros()
L1$imputeAll()

## Anova.R
L1$anovaAnalysisPreImpute()
L1$anovaAnalysis()
L1$drawVenDiagramSignificant()
L1$drawVolcanoPlot() #or run the following commands at the end (it might crash)


save(L1,file = "../backup.RData")


