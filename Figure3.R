setwd('path')

library(Seurat)
library(ggplot2)
library(data.table)


#Figure 3a
dat = readRDS('data/seurat_covid_pbmc_anno_upload.rds')
Idents(dat) = 'celltype'
NK_data = subset(dat,celltype %in% "NK cells")
dim(NK_data)

markers = c('CD3D','CD3G','CD3E','CD4','CD40LG','CD8A','CD8B',
            'GNLY','NKG7','NCAM1',"FCGR3A",
            "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DPA1","HLA-DPB1",
            'SELL','KLRB1','KLRC1',"KLRC2",'GZMA','GZMB','GZMH','GZMK','B3GAT1'
             )

NK_data$cell_cluster = factor(NK_data$cell_cluster,levels = rev(sort(unique(NK_data$cell_cluster))))
Idents(NK_data)='cell_cluster'


pdf('plot/Fig3a.pdf',height = 4,width = 12)
DotPlot(NK_data,features = markers, scale = T,assay = 'RNA')+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1))+
  theme(axis.text.x = element_text(angle = 90,hjust =1,vjust = 0.5,face = 'italic'))+
  labs(x='',y='')
dev.off()


#Figure 3b
metadata = dat@meta.data
data = subset(dat,cell_cluster %in% c("c37_NK-CD3-GZMK-TCF7", "c38_NK-CD3-GZMB-ZBTB16", 
                                      "c39_NK-CD3-GZMB-CX3CR1", "c40_NK-CD3-GZMB-ZNF683", "c41_NK-CD3-GZMB-LAG3", 
                                      "c42_NK-CD3-MKI67") | celltype %in% "T cells" )

#
exp <- FetchData(data,vars = c("CD3D","CD3E","CD3G"),slot = 'data') %>% as.data.frame()
exp <- as.data.frame(apply(exp, 2, function(x) ifelse(x > 0, "+", "-")))
exp <- merge(exp,metadata,by = 'row.names',all.x=T)
exp$cell_category = ifelse(exp$celltype == 'T cells','T cells','CD3+ NK cells')


stat <- exp %>% dplyr::group_by(cell_category,CD3D, CD3E, CD3G) %>% dplyr::summarise(count = n())
stat$proportion = proportions(stat$count)
stat$group = paste0(stat$CD3D,stat$CD3E,stat$CD3G)

# barplot
stat$group = factor(stat$group,levels = c("+++", "++-", "+-+", "-++", "+--","-+-", "--+", "---"))
color_used = c("#EB6F5D","#4DBBD5","#5050FF","#F39B7F","#66A61E","#F0E685","#EE4C97","#6376A0")

p=ggplot(stat) + 
  geom_bar(aes(x = cell_category,y = proportion,fill = group),
           stat = "identity",position = 'fill',width = 0.7,size = 0.5)+ 
  theme_classic() +
  coord_flip()+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = color_used)+
  labs(x='',y = 'Proportion',fill = '')+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        strip.background = element_rect(fill = "grey90", linetype = "solid"),
        axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 12),
        legend.position = "top"
  )+
  guides(color = guide_legend(keywidth=0.5, keyheight=0.5, default.unit="inch"))+
  theme(legend.key.size = unit(0.3, "inch"))

ggsave('plot/Fig3b.pdf',p,width = 4,height = 2.5)



#Figure 3c
metadata = dat@meta.data
TCR_data = read.table('data/TCR_data.tsv',sep = '\t',header = T)
TCR_data = merge(TCR_data,metadata,by = 'row.names')
TCR_data = subset(TCR_data,celltype %in% c("T cells","NK cells"))

TCR_data <- TCR_data %>%
  mutate(cell_category = case_when(
    cell_cluster %in% c("c15_CD4_Tn-CCR7", "c16_CD4_Treg-FOXP3", "c17_CD4_Tfh like-CXCR5-GPR183", 
                        "c18_CD4_Th1 like-TBX21-GZMB", "c19_CD4_Th1 like-CXCR3-GZMK", 
                        "c20_CD4_Th2 like-CCR4-GATA3", "c21_CD4_Th17 like-CCR6") ~ 'CD4+ T',
    cell_cluster %in% c("c22_CD8_Tn-CCR7", "c23_CD8_Tcm-GPR183", "c24_CD8_Tem-GZMK", 
                        "c25_CD8-GNLY-GZMK", "c26_CD8-ZNF683-FOS", "c27_CD8-GNLY-GZMB-CX3CR1", 
                        "c28_CD8-GNLY-GZMB-ZNF683", "c29_CD8-GNLY-GZMH", "c30_CD8-KLRB1-TIGIT", 
                        "c31_CD8-MKI67") ~ 'CD8+ T',
    cell_cluster == "c34_gdT" ~ 'gdT',
    cell_cluster == "c35_MAIT" ~ 'MAIT',
    cell_cluster %in% c("c37_NK-CD3-GZMK-TCF7", "c38_NK-CD3-GZMB-ZBTB16", 
                        "c39_NK-CD3-GZMB-CX3CR1", "c40_NK-CD3-GZMB-ZNF683", "c41_NK-CD3-GZMB-LAG3", 
                        "c42_NK-CD3-MKI67") ~ 'CD3+ NK',
    cell_cluster %in% c("c43_NK-FCGR3A", "c44_NK-NFKB", "c45_NK-GZMK", "c46_NK-MKI67") ~ 'Classical NK',
    TRUE ~ 'Others' 
  ))


TCR_data$TCR = ifelse(TCR_data$Frequency==0,'TCRαβ-','TCRαβ+')
plot_data <- TCR_data %>% group_by(celltype,cell_category,TCR) %>% dplyr :: summarise(count = n())
plot_data <- subset(plot_data,cell_category != 'Others' )
plot_data$cell_category = factor(plot_data$cell_category,levels = c("CD4+ T", "CD8+ T","gdT", "MAIT", "Classical NK", "CD3+ NK"))
plot_data$TCR = factor(plot_data$TCR,levels = c('TCRαβ+','TCRαβ-'))
plot_data$celltype = factor(plot_data$celltype,levels = c("T cells","NK cells" ))

p=ggplot(plot_data, aes(x=cell_category, y= count))+
  geom_bar(aes(fill = TCR),position = "fill",stat = 'identity',width = 0.8) +
  scale_y_continuous(labels = scales::percent)+
  facet_grid(.~celltype,scales = 'free_x',space = 'free_x')+
  labs(fill='',y='Proportion',x='')+
  scale_fill_manual(values = c('#6565FF','grey'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 8),
        axis.text = element_text(color = 'black'),
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.grid = element_blank(),
        #strip.text = element_text(angle = 45)
  )
ggsave('plot/Fig3c.pdf',p,width = 4,height = 2.8)



#Figure 3d
NK_df = NK_data@meta.data
NK_df = NK_df %>%
  mutate(cell_category = case_when(
    cell_cluster %in% c("c37_NK-CD3-GZMK-TCF7", "c38_NK-CD3-GZMB-ZBTB16", 
                        "c39_NK-CD3-GZMB-CX3CR1", "c40_NK-CD3-GZMB-ZNF683", "c41_NK-CD3-GZMB-LAG3", 
                        "c42_NK-CD3-MKI67") ~ 'CD3+ NK cells',
    cell_cluster %in% c("c43_NK-FCGR3A", "c44_NK-NFKB", "c45_NK-GZMK", "c46_NK-MKI67") ~ 'Classical NK cells',
    TRUE ~ 'Others' ))

df = NK_df %>% group_by(sample,cell_category) %>% dplyr::summarise(count =  n())
df <- df %>% group_by(sample)%>% 
  dplyr::mutate(total_count = sum(count),
                proportion = round(count / total_count, digits = 2))%>% 
  data.frame()

df$sample = factor(df$sample,levels = c(c("HD_279","HD_281","HD_324","V_A_53","V_P_34","V_P_68","V_S_54","V_S_56","VII_A_18","VII_P_16","VII_P_48","VII_S_11","VII_S_44")))
df$cell_category = factor(df$cell_category ,levels = c("Classical NK cells","CD3+ NK cells"))

p=ggplot(df) + 
  geom_bar(aes(x = sample,y = proportion,fill = cell_category),
           stat = "identity",position = 'fill',width = 0.7,size = 0.5)+ 
  geom_text(aes(x = sample, y = proportion, label=scales::percent(proportion)),
            stat = "identity", size = 3, color = "white") + #
  theme_classic() +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(name='NK cells cluster',values = c("#E16CB7","#A032CB" ))+
  labs(x='Sample',y = 'Proportion')+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        strip.background = element_rect(fill = "grey90", linetype = "solid"),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

ggsave(str_glue('plot/Fig3d.pdf'),p,width = 8,height = 3)
