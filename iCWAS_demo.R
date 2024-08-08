
setwd('path')


# install.packages("devtools")
# devtools::install_github("ZxZhou4150/Redeconve", build_vignettes = F)


rm(list=ls())
library(Redeconve)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)

source('iCWAS_function.R')

#Input
####input1:bulk RNA-seq data, row:gene,col:sample
bulk_data <- read.table('data/bulk_tpm.tsv',sep = '\t',header = T,row.names = 1) 

####input2:scRNA-seq cell cluster or single cell expression
sc_data <- read.table('data/sc_tpm.tsv',sep = '\t',header = T,row.names = 1,check.names = F)
sc_data <- sc_data %>% Matrix::as.matrix()

#input3:phenotype data
phenotype = read.table('data/phenotype.tsv',sep = '\t',header = T,check.names = F)

#ICWAS analysis
iCWAS_res = iCWAS(sctype = 'tpm',deconve_hp = 1e+09,variables = 'continuous',
                 sample_ID = 'sample',sample_group = NULL )
