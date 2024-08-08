setwd('path')

library(Seurat)
library(ggplot2)
library(data.table)


#Figure 2a
dat = readRDS('sc_seurat.rds')
Idents(dat) = 'celltype'
NK_data = subset(dat,celltype %in% "NK cells")
dim(NK_data)

#Figure S2a
NK_data$cell_cluster = factor(NK_data$cell_cluster,levels = rev(sort(unique(NK_data$cell_cluster))))
Idents(NK_data)='cell_cluster'

NK_data$label = tstrsplit(NK_data$cell_cluster, "_")[[1]] 
NK_data$label = gsub("c", "", NK_data$label)

pdf(str_glue('plot/FigS2a.pdf'),height = 4,width = 8)
DimPlot(NK_data,cols= c("#A032CB","#D59CBD","#FC8D62","#CCECFF","#FED9A6",
                        "#6276B2","#E16CB7","#908EBC","#C5D9DF","#DBEBF7"),
        raster = F,group.by = 'cell_cluster')+
  labs(title = '')+NoAxes()+
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()


pdf(str_glue('plot/FigS2a_label.pdf'),height = 4,width = 4.5)
DimPlot(NK_data,label = T,cols= c("#A032CB","#D59CBD","#FC8D62","#CCECFF","#FED9A6",
                                  "#6276B2","#E16CB7","#908EBC","#C5D9DF","#DBEBF7"),
        raster = F,group.by = 'label',
        label.size = 3,repel = T)+labs(title = '')+
  NoAxes()+theme(legend.position = 'none')
dev.off()



#Figure S2b
markers = c('CD3D','CD3G','CD3E','CD4','CD40LG','CD8A','CD8B',
            'GNLY','NKG7','NCAM1',"FCGR3A",
            "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DPA1","HLA-DPB1",
            'SELL','KLRB1','KLRC1',"KLRC2",'GZMA','GZMB','GZMH','GZMK','B3GAT1'
)

pdf('plot/FigS2b.pdf',height = 13,width = 14)
FeaturePlot(NK_data, features = markers,ncol = 5)
dev.off()

