
setwd('path')

library(Seurat)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)


#Figure S6a
dat = readRDS('sc_seurat.rds')
Idents(dat) = 'celltype'
NK_data = subset(dat,celltype %in% "NK cells")

NK_data$cell_category =  case_when(
  NK_data$cell_cluster %in% c("c37_NK-CD3-GZMK-TCF7", "c38_NK-CD3-GZMB-ZBTB16", 
                              "c39_NK-CD3-GZMB-CX3CR1", "c40_NK-CD3-GZMB-ZNF683", "c41_NK-CD3-GZMB-LAG3", 
                              "c42_NK-CD3-MKI67") ~ 'CD3+ NK cells',
  NK_data$cell_cluster %in% c("c43_NK-FCGR3A", "c44_NK-NFKB", "c45_NK-GZMK", "c46_NK-MKI67") ~ 'Classical NK cells' )

NK_data$cell_category = factor(NK_data$cell_category ,levels = c("Classical NK cells","CD3+ NK cells"))


#Figure S6a
Idents(NK_data) = 'cell_category'
markers = FindMarkers(NK_data,ident.1 = "CD3+ NK cells", ident.2 = 'Classical NK cells',
                      logfc.threshold = 0.25,  min.pct = 0.1)
markers$gene = row.names(markers)
markers = dplyr::select(markers,gene,everything())
write.table(markers,'data/NK_cell_category_marker.tsv',sep = '\t',row.names = F)


#KEGG enrichment
gene_all <- bitr(rownames(dat), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db,drop = F)

CD3NK_up = subset(markers,avg_log2FC > 0)
CD3NK_up_ENTREZID = bitr(CD3NK_up$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db,drop = F) 

KEGG <- enrichKEGG(gene = CD3NK_up_ENTREZID$ENTREZID,
                   organism='hsa',universe=gene_all$ENTREZID,
                   pvalueCutoff=0.05, qvalueCutoff=0.05) 
KEGG = setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df = KEGG@result

df = df[order(df$p.adjust),] %>% head(10)
df$Description = factor(df$Description,levels = rev(df$Description))

p=ggplot(df, aes(x = Description, y = -log10(p.adjust))) +
  geom_bar(aes(fill=-log10(p.adjust)),stat = 'identity')+
  theme_classic()+
  coord_flip()+
  scale_fill_continuous(low="#F6A392",high="#7F0000")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
  labs(x='')+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(size = 12,margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(size = 12,margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5,
        legend.position = 'none')

ggsave('plot/FigS6.pdf',p,width = 5,height = 4)
  
