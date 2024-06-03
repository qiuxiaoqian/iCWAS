setwd('/gold160t/qiuxiaoqian/SARS_PBMC/Figures/github_upload')


# install.packages("devtools")
# devtools::install_github("ZxZhou4150/Redeconve", build_vignettes = F)


rm(list=ls())
library(Redeconve)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)



#Input
####input1:bulk RNA-seq data，row:gene,col:sample
conve_data <- read.table('data/bulk_TPM_table.tsv',sep = '\t',header = T,row.names = 1) 
                         

####input2:scRNA-seq cell cluster or single cell expression
sref <- read.table('data/sc_cell_cluster_readscount.tsv',sep = '\t',header = T,row.names = 1,check.names = F)
sref <- sref  %>% Matrix::as.matrix()


#input3:phenotype data
ELISPLOT_res = openxlsx::read.xlsx('data/ELISPOT.xlsx',rowNames = T,check.names = F)
ELISPLOT_res = ELISPLOT_res[names(conve_data),]

log_ELISPLOT_res = log2(ELISPLOT_res+1)



#Deconvolution
gene.list = intersect(rownames(conve_data),rownames(sref))
print(str_glue("get filtered genes:",length(gene.list)))

cell_abundance <- deconvoluting(sref, conve_data, genemode= "default",gene.list = gene.list, 
                     hpmode = "customized",hp = 1e+09,dopar = T, ncores = 64,realtime = F)

cell_abundance <- as.matrix(cell_abundance)
#cell_abundance <- cell_abundance[order(rownames(cell_abundance)),sample]

cell_abundance <- apply(cell_abundance, 2, function(x) x / sum(x) *100) #Proportion
write.table(cell_abundance,'data/bulk_estimated_cell_proprotion.tsv',sep = '\t')



#Association analysis
library(psych)
cell_abundance = t(cell_abundance) %>% data.frame(check.names = F)
pearson <- corr.test(cell_abundance,log_ELISPLOT_res,use="pairwise", method="pearson", 
                     adjust="BH", alpha=.05, ci=TRUE, minlength=10)


#import result
pearson_r <- pearson$r %>% melt()
pearson_r  <- unite (pearson_r,"pair","Var1","Var2", sep=":", remove = T) 
pearson_p <- pearson$p %>% melt()
pearson_p  <- unite (pearson_p,"pair","Var1","Var2", sep=":", remove = T) 
pearson_padj <- pearson$p.adj %>% melt()
pearson_padj  <- unite (pearson_padj,"pair","Var1","Var2", sep=":", remove = T) 

pearson_cor <- left_join(pearson_r,pearson_p,by="pair")
pearson_cor <- left_join(pearson_cor,pearson_padj,by="pair")
names(pearson_cor) <- c('pair','r','p','padj')

write.table(pearson_cor,"data/iCWAS_result.tsv",sep = '\t',row.names = F)

