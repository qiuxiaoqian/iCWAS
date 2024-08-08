setwd('path')


library(Seurat)
library(ggplot2)
library(data.table)


#Figure 1b
dat = readRDS('path/sc_seurat.rds')

color_used = c("#A6CEE3","#1F78B4","#FB9A99","#B2DF8A","#E69C8D","#F07590","#33A02C","#FDBF6F",
               "#66C2A5","#F0E685","#E31A1C","#8DA0CB","#E78AC3","#ECB888","#B589BC","#91CDC1",
               "#77A1EC","#EAD7EC","#A3C4A5","#EED785","#CD7560","#EFE0E7","#7FC28E","#EEBED6",
               "#F7D39B","#EDB869","#D4B595","#E3EDC0","#F7E3DB","#F48F89","#F9C9BD","#84D7F7",
               "#D4EEF6","#FEF6AE","#C76DA8","#FC7440","#A032CB","#D59CBD","#FC8D62","#CCECFF",
               "#FED9A6","#6276B2","#E16CB7","#908EBC","#C5D9DF","#DBEBF7","#85C1D7","#AF88BB",
               "#DFB9D5","#B099B5","#DEDBEE","#5394C3","#B8A89F","#917393","#F5BCC9","#5DB4F3",
               "#F7F4B7","#6F89D3","#EA94CD","#D8B6F7","#CCEBC5","#DECBE4","#E5D8BD","#F4CAE4",
               "#DCDCDC","#B3E2CD","#CBD5E8","#E6F5C9","#FFF2AE","#BEAED4")

#legend
pdf('plot/Fig1b.pdf',height = 7,width = 15)
DimPlot(dat,cols = color_used,raster = F,group.by = 'cell_cluster')+
  labs(title = '')+NoAxes()+
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 3))
dev.off()


#label
options(ggrepel.max.overlaps = Inf)
dat$label = tstrsplit(dat$cell_cluster, "_")[[1]] 
dat$label = gsub("c", "", dat$label)

pdf('plot/Fig1b_label.pdf',height = 5.5,width = 5.5)
DimPlot(dat,label = T,cols = color_used,raster = F,group.by = 'label',
        label.size = 3,repel = T)+labs(title = '')+
  NoAxes()+theme(legend.position = 'none')
dev.off()




#Figure 1c
library(ComplexHeatmap)
library(RColorBrewer)

df = read.table('data/bulk_tpm.tsv',sep = '\t',header = T)
metadata = read.csv('data/bulk_sample_metadata.csv')
exp.df = t(df) %>% data.frame()
exp.df = merge(metadata,exp.df,by.x = 'sample',by.y = 'row.names')


#ANOVA to select plot gene:p_value < 0.00001 $ cell cluster marker gene
calculate_annova_more <- function(data_frame,factor_list){
  c_fomula <- paste(gene, " ~ ",paste(factor_list,collapse = '*'),sep="")
  c_table <- summary(aov(as.formula(c_fomula),data=data_frame))
  c_table <- as.matrix(c_table[[1]])
  p_value_list <- c_table[1:(nrow(c_table)-1),5]
  p_value_list <- data.frame(t(p_value_list))
  names(p_value_list) = 'p_value'
  rownames(p_value_list)= gene
  return(p_value_list)
}

res_df = data.frame()
for (i in 4:ncol(exp.df)){
  gene = names(exp.df)[i]
  tryCatch({
    res = calculate_annova_more(data_frame = exp.df,factor_list = c("group"))
    res_df = rbind(res_df,res)
  }, error = function(e) {
    message(paste("Error occurred:", conditionMessage(e)))
    print('error and next')
  })
}

res_df$gene = row.names(res_df)
res_df = select(res_df,'gene',everything())

marker = c("XBP1","MS4A1","CD19","FCGR3A","NCAM1",
           "CD3D","CD3E","CD3G","CD4","CD40LG","CD8A","CD8B","TRDV2","TRGV9","SLC4A10",
           "CST3","LYZ","CLEC9A","CD1C","LILRA4","CD68","KIT","PPBP","MKI67")

res_df = subset(res_df,p_value < 0.00001 | gene %in% marker)

GOI = intersect(row.names(df),res_df$gene)
GOI_df = df[GOI,]
df_scaled <- t(scale(t(GOI_df))) %>% data.frame()
#write.csv(df_scaled,'plot/plot_data/Fig1c.csv')

#heatmap plot
index <- which(rownames(df_scaled) %in% marker)
labs <- rownames(df_scaled)[index]

color_mapping = setNames(c("#82A5D6","#51A246","#6E51A0"), c("Healthy", "BA.5", "BF.7"))
lab1 = columnAnnotation(Groups = c(rep("Healthy",9),rep("BA.5",15),rep("BF.7",14)),
                        show_annotation_name = F,
                        annotation_name_gp = gpar(fontsize = 4),
                        col = list(Groups = color_mapping))

lab2 = rowAnnotation(foo = anno_mark(at = index,
                                     labels = labs,
                                     labels_gp = gpar(fontsize = 12),
                                     lines_gp = gpar(lwd=0.2)))



col_fun = colorRampPalette(brewer.pal(9, "Oranges"))(50)

pdf("plot/Fig1c.pdf",width = 10, height = 6.2)
Heatmap(df_scaled,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        col = col_fun,
        column_title = "Samples",column_title_side = "bottom", #标题位置
        column_title_gp = gpar(fontsize = 20,col = "black"),
        cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE,show_column_dend = FALSE,
        show_row_names = FALSE,
        bottom_annotation = lab1,right_annotation = lab2,
        name = "Normalized expression")

dev.off()



#Figure 1d
cell_abundance = read.table('data/bulk_estimated_cell_proprotion.tsv',sep = '\t')
cell_df <- as.matrix(cell_abundance)
cell_df <- melt(cell_df)
colnames(cell_df) <- c("cell_cluster","sample","proportion")


#
celltype_metadata <- read.table('data/sc_celltype_metadata.tsv',sep = '\t',header = T)
cell_df <- left_join(cell_df,celltype_metadata,by="cell_cluster")

color_used = c("#A6CEE3","#E69C8D","#B589BC","#EFE0E7","#84D7F7",
               "#A032CB","#85C1D7","#EA94CD","#E5D8BD","#E6F5C9")
cell_df$sub_celltype= factor(cell_df$sub_celltype,levels = c("Plasma cells", "B cells","CD4+ T cells", "CD8+ T cells","Other T cells",
                                                              "NK cells","Monocytes","Dendritic cells","Other myeloid cells","Platelet" ) )

metadata = read.csv('data/bulk_sample_metadata.csv')
cell_df = left_join(cell_df,metadata,by='sample')
cell_df$sample = factor(cell_df$sample,levels = unique(cell_df$sample))

p=ggplot(cell_df) + 
  geom_bar(aes(x = sample,y = proportion,fill = sub_celltype),
           stat = "identity",position = 'stack',width = 0.7,size = 0.5)+ 
  theme_classic() +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = color_used)+
  labs(x='Samples',y = 'Proportion',fill = NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.2, linetype="solid"),
        strip.background = element_rect(fill = "grey90", linetype = "solid"),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
  )+
  theme(legend.key.size = unit(0.25, "inch"))

ggsave('plot/Fig1d.pdf',p,width = 7,height = 4)




#Figure 1e

#iCWAS analysis


#Input
####input1:bulk RNA-seq data, row:gene,col:sample
bulk_data <- read.table('data/bulk_tpm.tsv',sep = '\t',header = T,row.names = 1) 

####input2:scRNA-seq cell cluster or single cell expression
sc_data <- read.table('data/sc_cell_cluster_readscount.tsv',sep = '\t',header = T,row.names = 1,check.names = F)
sc_data <- sc_data %>% Matrix::as.matrix()

#input3:phenotype data
phenotype = openxlsx::read.xlsx('data/ELISPOT.xlsx',rowNames = T,check.names = F)
phenotype = phenotype[names(bulk_data),]
phenotype = log2(phenotype+1)



#ICWAS analysis
for (variant in names(phenotype)[2:ncol(phenotype)]){
  
  sub_phenotype = phenotype[,c(sample,variant)]
  
  iCWAS_res = iCWAS(bulk_data,sc_data,sub_phenotype,
                    sctype = 'tpm',deconve_hp = 1e+09,variables = 'continuous',
                    sample_ID = 'sample',sample_group = NULL )
  
  Scissor(bulk_dataset=bulk_data,sc_dataset=sc_data,sub_phenotype,
          tag = tag, alpha = 0.4, save_file = str_glue("{variant}_Scissor_res"),family = "binomial")
  
  
}
  


#plot for iCWAS or Scissor
metadata = dat@meta.data


#Add matadata and Plot


#iCWAS
data = read.table('data/iCWAS_result.tsv',sep = '\t',header = T)
data = separate(data,col='pair',into = c('cell_cluster','variant'),sep = ':')
data = subset(data,variant %in% c("BA.5/BF.7","SARS-COV-2"))
data$label = tstrsplit(data$cell_cluster, "_")[[1]] 


celltype_metadata <- read.table('data/sc_celltype_metadata.tsv',sep = '\t',header = T)
data <- left_join(data,celltype_metadata,by=c('cell_cluster'))

data$variant = factor(data$variant,levels = c("SARS-COV-2","BA.5/BF.7"))
data$sub_celltype= factor(data$sub_celltype,levels = c("Plasma cells", "B cells","CD4+ T cells", "CD8+ T cells","Other T cells",
                                                       "NK cells","Monocytes","Dendritic cells","Other myeloid cells","Platelet" ) )

p = ggplot(data, aes(x = label, y = variant)) +
  geom_point(aes(size=abs(r),color=r),shape=16,alpha = 0.9)+#
  labs(x="",  y="")+
  theme_bw()+
  scale_size_continuous(range = c(1, 5))+
  scale_color_gradient2(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF",#low = "#424da7", mid = "#ffffff", high = "#dd2b19",
                        breaks = c(-0.4, -0.2, 0, 0.2, 0.4), 
                        labels = c("-0.4", "-0.2", "0", "0.2", "0.4"))+
  theme(axis.text = element_text(size = 8, color="black"),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="black", linewidth =0.5, linetype="solid"),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 10)
  )


col= c("#A6CEE3","#E69C8D","#B589BC","#EFE0E7","#84D7F7",
       "#A032CB","#85C1D7","#EA94CD","#E5D8BD","#E6F5C9")


theme_niwot <- function(){
  theme(
    legend.key=element_blank(), 
    legend.text = element_text(color="black",size=10),
    legend.spacing.x=unit(0.1,'cm'),
    legend.key.width=unit(0.5,'cm'),
    legend.key.height=unit(0.5,'cm'), 
    legend.background=element_blank())
}


p_anno=data %>% select(6,8) %>% mutate(group="sub_celltype") %>% 
  as_tibble() %>% 
  mutate(sub_celltype=as.character(sub_celltype)) %>% 
  ggplot(aes(label,group,fill=sub_celltype))+
  geom_tile()+
  scale_fill_manual(values=col)+
  scale_y_discrete(expand = c(0,0),position="right")+
  scale_x_discrete(expand=c(0,0))+
  theme_void()+
  #theme(axis.text.y=element_text(color="black",size=10))+
  theme_niwot()+
  theme(legend.position = 'none')


library(aplot)
p=p %>% insert_bottom(p_anno,height = 0.005)
ggsave('Fig1e_iCWAS.pdf',p,width = 7.945,height = 20)




#Scissor
infos = readRDS(str_glue("{variant}_Scissor_res"))
Scissor_select <- rep('Background', ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos$Scissor_pos] <- 'Scissor+'
Scissor_select[infos$Scissor_neg] <- 'Scissor-'
Scissor_select <- data.frame(Scissor_select)
Scissor_select$celltype = rownames(Scissor_select)
names(Scissor_select) = c(pro, "celltype")


meta = left_join(metadata,Scissor_select,by ='cell_ids')
da <- meta %>%
  group_by(cell_cluster, Variant, value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(cell_cluster, Variant, value, fill = list(count = 0)) %>%
  group_by( cell_cluster, Variant) %>%
  mutate(total_count = sum(count),
         percentage = count / total_count) %>%
  ungroup() %>%
  data.frame()

da = subset(da,value %in% c('Scissor+','Scissor-'))
p = ggplot(da, aes(x = label, y = value)) +
  geom_point(aes(size=abs(percentage),color=percentage),shape=16,alpha = 0.9)+#
  labs(x="",  y="")+
  theme_bw()+
  scale_size_continuous(range = c(0, 5))+
  scale_color_gradient2(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF",#low = "#424da7", mid = "#ffffff", high = "#dd2b19",
                        breaks = c(0, 0.25, 0.5, 0.5, 0.75, 1),
                        labels = c(0, 0.25, 0.5, 0.5, 0.75, 1))+
  theme(axis.text = element_text(size = 8, color="black"),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="black", linewidth =0.5, linetype="solid"),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 10),
        strip.background = element_blank(),
        #strip.text = element_blank()
  )


p_anno= ggplot(da,aes(x=label,y=0.5,fill = cell_class4))+
  geom_tile()+
  scale_fill_manual(values=c("#A6CEE3","#E69C8D","#B589BC","#EFE0E7","#84D7F7",
                             "#A032CB","#85C1D7","#EA94CD","#E5D8BD","#E6F5C9"))+
  scale_y_discrete(expand = c(0,0),position="right")+
  scale_x_discrete(expand=c(0,0))+
  theme_void()+
  #theme(axis.text.y=element_text(color="black",size=10))+
  theme(
    legend.key=element_blank(),  
    legend.text = element_text(color="black",size=10), 
    legend.spacing.x=unit(0.1,'cm'), 
    legend.key.width=unit(0.5,'cm'), 
    legend.key.height=unit(0.5,'cm'), 
    legend.background=element_blank(),
    legend.position = 'none')


library(aplot)
p=p %>% insert_bottom(p_anno,height = 0.005)
ggsave('Fig1e_Scissor.pdf',p,width = 7.945,height = 20)

