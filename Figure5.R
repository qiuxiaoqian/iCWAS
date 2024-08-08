setwd('path')

library(Seurat)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)


#Figure 5a
dat = readRDS('sc_seurat.rds')
Idents(dat) = 'celltype'
NK_data = subset(dat,celltype %in% "NK cells")

NK_data$cell_category =  case_when(
  NK_data$cell_cluster %in% c("c37_NK-CD3-GZMK-TCF7", "c38_NK-CD3-GZMB-ZBTB16", 
                        "c39_NK-CD3-GZMB-CX3CR1", "c40_NK-CD3-GZMB-ZNF683", "c41_NK-CD3-GZMB-LAG3", 
                        "c42_NK-CD3-MKI67") ~ 'CD3+ NK cells',
  NK_data$cell_cluster %in% c("c43_NK-FCGR3A", "c44_NK-NFKB", "c45_NK-GZMK", "c46_NK-MKI67") ~ 'Classical NK cells' )

NK_data$cell_category = factor(NK_data$cell_category ,levels = c("Classical NK cells","CD3+ NK cells"))


pdf(str_glue('plot/Fig5a-1.pdf'),height = 4,width = 4)
DimPlot(NK_data,cols= c("#E16CB7","#A032CB" ),
        raster = F,group.by = 'cell_category')+
  labs(title = '')+NoAxes()+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.position = "bottom")
dev.off()


pdf('plot/Fig5a-2.pdf',height = 6.5,width = 5.6)
FeaturePlot(NK_data, features = c("HLA-DRA","HLA-DRB1","HLA-DRB5",
                                  "HLA-DQA1","HLA-DPA1","HLA-DPB1"),ncol = 2)
dev.off()


#Figure 5b
Idents(NK_data) = 'cell_category'
markers = FindMarkers(NK_data,ident.1 = "CD3+ NK cells", ident.2 = 'Classical NK cells',
                      logfc.threshold = 0.25,  min.pct = 0.1)
markers$gene = row.names(markers)
markers = dplyr::select(markers,gene,everything())
write.table(markers,'data/NK_cell_category_marker.tsv',sep = '\t',row.names = F)


#GO enrichment
CD3NK_up = subset(markers,avg_log2FC > 0)
CD3NK_up_ENTREZID = bitr(CD3NK_up$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db,drop = F) 

GO <- enrichGO(gene = CD3NK_up_ENTREZID$ENTREZID,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff = 0.05, qvalueCutoff = 0.05,readable = T) 

GO = setReadable(GO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df = GO@result

df = df[order(df$p.adjust),] %>% head(10)
df$Description = factor(df$Description,levels = rev(df$Description))

p=ggplot(df, aes(x = Description, y = -log10(p.adjust))) +
  geom_bar(aes(fill=-log10(p.adjust)),stat = 'identity',width = 0.7)+
  theme_classic()+
  coord_flip()+
  scale_fill_continuous(low="#F6A392",high="#7F0000")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
  labs(x='')+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(size = 12,margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(size = 12,margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        legend.title.align=0.5,
        legend.position = 'none')

ggsave('plot/Fig5b.pdf',p,width =5,height = 4)
  