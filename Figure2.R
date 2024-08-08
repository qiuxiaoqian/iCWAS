

#Profiling cell subpopulations associated with tumor immunotherapy and autoimmune disease
setwd('path/Melanoma_data')

rm(list=ls())
library(Redeconve)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)

devtools::install_github('sunduanchen/Scissor')
library(Scissor)

source('iCWAS_function.R')


#Melanoma Dadaset

#Input
####input1:bulk RNA-seq data, row:gene,col:sample
bulk_data <- read.table('bulk_tpm.tsv',sep = '\t',header = T,row.names = 1) 

####input2:scRNA-seq cell cluster or single cell expression
sc_data <- read.table('sc_tpm.tsv',sep = '\t',header = T,row.names = 1,check.names = F)
sc_data <- sc_data %>% Matrix::as.matrix()

#input3:phenotype data
phenotype = read.table('phenotype.tsv',sep = '\t',header = T,check.names = F)


#ICWAS analysis
iCWAS_res = iCWAS(bulk_data,sc_data,phenotype,sctype = 'tpm',deconve_hp = 1e+6,variables = 'discrete',
                  sample_ID = "Patient_ID",sample_group = "Response_Status")



#Scissor analysis
tag =c("Non-Responding","Responding")
phenotype[,"Response_Status"] = plyr::mapvalues(phenotype[,"Response_Status"],tag,c(0,1))

Scissor_res = Scissor(bulk_dataset=bulk_data,sc_dataset=sc_data, phenotype ,
        tag = tag, alpha = 0.05, save_file = "Scissor_res",family = "binomial")
        



#Add matadata and Plot
dat = readRDS('seurat_obj.rds')
metadata = dat@meta.data

#Add Scissor result
Scissor_res = readRDS('Scissor_res.rds')
Scissor_select <- rep('Background', ncol(sc_data))
names(Scissor_select) <- colnames(sc_data)
Scissor_select[Scissor_res$Scissor_neg] <- 'Scissor-'
Scissor_select[Scissor_res$Scissor_pos] <- 'Scissor+'
Scissor_select <- data.frame(Scissor_select)
Scissor_select$cell_ids = rownames(Scissor_select)
meta = left_join(metadata,Scissor_select,by ='cell_ids')


#Add iCWAS result
estimated_cell_pro = read.table(str_glue('bulk_estimated_cell_proprotion.tsv'),header = T)
estimated_cell_pro = estimated_cell_pro[rowSums(estimated_cell_pro != 0) > 0, ]
estimated_cell_pro = tibble::rownames_to_column(estimated_cell_pro,'cell_ids')

iCWAS_result = read.table('iCWAS_result.tsv',header = T)
names(iCWAS_result)[4] = 'log2FC'

iCWAS_result = subset(iCWAS_result,cell %in% estimated_cell_pro$cell_ids)

data = iCWAS_result %>% 
  mutate(sig = case_when(log2FC>0 & p.value<0.05 ~ "Up",log2FC<0 & p.value< 0.05 ~ "Down",TRUE ~ "NotSig")) %>% 
  mutate(iCWAS = case_when(log2FC > 0 ~ "iCWAS+",log2FC < 0 ~ "iCWAS-",log2FC == 0  ~ "Background"))

meta = left_join(meta,unique(data[,c("cell","iCWAS")]),by = c('cell_ids'='cell'))
meta$iCWAS = ifelse(is.na(meta$iCWAS),'Background',meta$iCWAS)

rownames(meta) = meta$cell_ids

dat@meta.data = meta
dat$Scissor_select = factor(dat$Scissor_select,levels = c("Scissor+","Scissor-","Background" ))
dat$iCWAS = factor(dat$iCWAS,levels = c( "iCWAS+" ,"iCWAS-","Background"))


dat$NonMalignantCellType = factor(dat$NonMalignantCellType,
                                  levels = c("B", "T", "NK", "Malignant", "CAF", "Macro", "Endo"))

pdf(str_glue('dat_CellType_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "NonMalignantCellType",label = F)+
  labs(title = str_glue('Cell type'),x='UMAP1',y='UMAP2')
dev.off()

pdf(str_glue('dat_Scissor_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "Scissor_select",label = F,
        cols = c("Scissor+"="red","Scissor-"="blue","Background"="#D9D9D9"))+
  labs(title = str_glue('Scissor'),x='UMAP1',y='UMAP2')
dev.off()

pdf(str_glue('dat_iCWAS_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "iCWAS",label = F,
        cols = c("iCWAS+"="red","iCWAS-"="blue","Background"="#D9D9D9"))+
  labs(title = str_glue('iCWAS'),x='UMAP1',y='UMAP2')
dev.off()





#SLE dataset

#Input
####input1:bulk RNA-seq data, row:gene,col:sample
bulk_data <- read.table('bulk_tpm.tsv',sep = '\t',header = T,row.names = 1) 

####input2:scRNA-seq cell cluster or single cell expression
sc_data <- read.table('sc_tpm.tsv',sep = '\t',header = T,row.names = 1,check.names = F)
sc_data <- sc_data %>% Matrix::as.matrix()

#input3:phenotype data
phenotype = read.table('phenotype.tsv',sep = '\t',header = T,check.names = F)


#ICWAS analysis
iCWAS_res = iCWAS(bulk_data,sc_data,phenotype,sctype = 'tpm',deconve_hp = 1e+6,variables = 'discrete',
                  sample_ID = "Sample",sample_group =  "Characteristics")


#Scissor analysis
tag = c("inactive","active")
phenotype[,"Characteristics"] = plyr::mapvalues(phenotype[,"Characteristics"],tag,c(0,1))

for (alpha in c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)){
  Scissor_res = Scissor(bulk_dataset=bulk_data,sc_dataset=sc_data, phenotype ,
                        tag = tag, alpha = 0.05, save_file = "alpha",family = "binomial")
  
}


#Add matadata and Plot
dat = readRDS('seurat_obj.rds')
metadata = dat@meta.data


#Scissor result
Scissor_res = readRDS('Scissor_res.rds')
Scissor_select <- rep('Background', ncol(sc_data))
names(Scissor_select) <- colnames(sc_data)
Scissor_select[Scissor_res$Scissor_neg] <- 'Scissor-'
Scissor_select[Scissor_res$Scissor_pos] <- 'Scissor+'
Scissor_select <- data.frame(Scissor_select)
Scissor_select$cell_ids = rownames(Scissor_select)
meta = left_join(metadata,Scissor_select,by ='cell_ids')


#iCWAS result
estimated_cell_pro = read.table(str_glue('bulk_estimated_cell_proprotion.tsv'),header = T)
estimated_cell_pro = estimated_cell_pro[rowSums(estimated_cell_pro != 0) > 0, ]
estimated_cell_pro = tibble::rownames_to_column(estimated_cell_pro,'cell_ids')

iCWAS_result = read.table('iCWAS_result.tsv',header = T)
names(iCWAS_result)[4] = 'log2FC'

iCWAS_result = subset(iCWAS_result,cell %in% estimated_cell_pro$cell_ids)

data = iCWAS_result %>% 
  mutate(sig = case_when(log2FC>0 & p.value<0.05 ~ "Up",log2FC<0 & p.value< 0.05 ~ "Down",TRUE ~ "NotSig")) %>% 
  mutate(iCWAS = case_when(log2FC > 0 ~ "iCWAS+",log2FC < 0 ~ "iCWAS-",log2FC == 0  ~ "Background"))


meta = left_join(meta,unique(data[,c("cell","iCWAS")]),by = c('cell_ids'='cell'))
meta$iCWAS = ifelse(is.na(meta$iCWAS),'Background',meta$iCWAS)

rownames(meta) = meta$cell_ids

dat@meta.data = meta
dat$Scissor_select = factor(dat$Scissor_select,levels = c("Scissor+","Scissor-","Background" ))
dat$iCWAS = factor(dat$iCWAS,levels = c( "iCWAS+" ,"iCWAS-","Background"))


dat$anno1 = factor(dat$anno1,levels = c("B", "CD4T","CD8T","MKI67_ProT", "NK", "CD14_Mono",  "CD16_Mono", 
                                        "LDG", "Mega", "cDC2", "pDC","Platelet")) 


color_used = c("#E64B35","#F39B7F","#66A61E","#FF7F00","#E6AB02","#4DBBD5",
               "#5050FF","#EE4C97","#7FC97F","#F0E685","#CCECFF","#D8B6F7")
pdf(str_glue('dat_CellType_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "anno1",label = F,cols = color_used)+
  labs(title = str_glue('Cell type'),x='UMAP1',y='UMAP2')
dev.off()

pdf(str_glue('dat_Scissor_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "Scissor_select",label = F,
        cols = c("Scissor+"="red","Scissor-"="blue","Background"="#D9D9D9"))+
  labs(title = str_glue('Scissor'),x='UMAP1',y='UMAP2')
dev.off()

pdf(str_glue('dat_iCWAS_umap.pdf'),height = 4,width = 5.5)
DimPlot(dat, reduction = "umap",group.by  = "iCWAS",label = F,
        cols = c("iCWAS+"="red","iCWAS-"="blue","Background"="#D9D9D9"))+
  labs(title = str_glue('iCWAS'),x='UMAP1',y='UMAP2')
dev.off()