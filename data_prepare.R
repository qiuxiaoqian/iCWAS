



#Figure 1b
metadata = dat@meta.data
write.table(metadata,'data/sc_cell_metadata.tsv',sep = '\t',row.names = F)


#bulk
sample = c("HD248", "HD258", "HD277", "HD279", "HD281", "HD283", "HD320","HD324", "HD247",
           "FA51", "FA53", "FP21", "FP34", "FP57", "FP67", "FP68", "FP69", 
           "FP71", "FS45", "FS48", "FS54", "FS56", "FS64", "FS70",  
           "SA18", "SA42", "SP16", "SP17", "SP20", "SP48", "SP53", 
           "SP64", "SS11", "SS24", "SS44", "SS47", "SS56", "SS60")

re_sample = c("HD_248", "HD_258", "HD_277", "HD_279", "HD_281", "HD_283", "HD_320","HD_324", "HD_247",
              "BTI_V51", "BTI_V53", "BTI_V21", "BTI_V34", "BTI_V57", "BTI_V67", "BTI_V68", "BTI_V69", 
              "BTI_V71", "BTI_V45", "BTI_V48", "BTI_V54", "BTI_V56", "BTI_V64", "BTI_V70",  
              "BTI_VII18", "BTI_VII42", "BTI_VII16", "BTI_VII17", "BTI_VII20", "BTI_VII48", "BTI_VII53", 
              "BTI_VII64", "BTI_VII11", "BTI_VII24", "BTI_VII44", "BTI_VII47", "BTI_VII56", "BTI_VII60")

df = read.table('/gold160t/qiuxiaoqian/SARS_PBMC/bulkRNA_seq/FeatureCounts/TPM_table.tsv',
                sep = '\t',header = T,row.names = 1) 

#重命名
names(df)
df = df[,sample]
names(df) = re_sample
names(df)

write.table(df,'data/bulk_TPM_table.tsv',sep = '\t')



df = read.table('data/bulk_TPM_table.tsv',sep = '\t',header = T)
df = tibble::column_to_rownames(df,var = 'Gene_name')
df = df[,sample]
names(df) = re_sample

#iCWAS pipeline
####input2:scRNA-seq cell cluster or single cell expression

sc = readRDS('/gold160t/qiuxiaoqian/SARS_PBMC/scRNA_seq/anno2/seurat_covid_pbmc_anno_upload.rds')
sc_list <- SplitObject(sc,split.by = 'cell_cluster')

sc_cell_cluster <- map(sc_list,function(x){
  counts <- x@assays$RNA@counts %>% Matrix::as.matrix()
  counts <- apply(counts, 1, sum)
  re <- data.frame(row.names = names(counts),exp = counts,stringsAsFactors = F)
  colnames(re) <- x$cell_cluster %>% unique
  return(re)
})


sc_cell_cluster <- sc_cell_cluster %>% do.call(cbind,.)
sc_cell_cluster$Gene_name <- rownames(sc_cell_cluster)
sc_cell_cluster <- dplyr::select(sc_cell_cluster,'Gene_name',everything())
#write.table(sc_cell_cluster,'data/sc_cell_cluster_count.tsv',sep = '\t',row.names = F)


metadata = read.table('data/sc_cell_metadata.tsv',sep = '\t',header = T)
celltype_metadata = metadata %>% select(.,c("celltype", "sub_celltype", "cell_cluster")) %>% distinct() %>% arrange(cell_cluster)
write.table(celltype_metadata,'data/sc_celltype_metadata.tsv',row.names = F,sep = '\t')



#添加样本信息
sample_metadata = read.csv(sample_metadata_file)



