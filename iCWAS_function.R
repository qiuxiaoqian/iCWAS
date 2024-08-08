


iCWAS = function(bulk_data,sc_data,phenotype,
                 sctype = 'count',deconve_hp = 1e+09,variables = 'continuous',
                 sample_ID = 'sample',sample_group = NULL ){
  
  #parameter use in iCWAS
  #' @bulk_data bulk RNA-seq expression matrix
  #' @sc_data single-cell expression matrix
  #' @phenotype an immune response vector paired with the bulk RNA-seq samples
  #' @sctype The quantification of scRNA-seq, either 'count' or 'TPM', gene*cell
  #' @deconve_hp The hyperparameter you want to use in deconvolution
  #' @variables The variable type of phenotype, either 'discrete' or 'continuous'
  #' @sample_ID The colname of sample in phenotype
  #' @sample_group The colname of sample group in phenotype
  
  
  
  #Step1:Deconvolution
  gene.list = intersect(rownames(bulk_data),rownames(sc_data))
  print(str_glue("get filtered genes:",length(gene.list)))
  
  if (sctype == 'count'){normalize = T}else{normalize = F}
  
  cell_abundance <- deconvoluting(sc_data, bulk_data, genemode= "default",gene.list = gene.list, 
                                  hpmode = "customized",hp = deconve_hp,normalize = normalize,
                                  dopar = T, ncores = 64,realtime = F)
  
  cell_abundance <- as.matrix(cell_abundance)
  
  cell_abundance <- apply(cell_abundance, 2, function(x) x / sum(x) *100) #Proportion
  write.table(cell_abundance,'bulk_estimated_cell_proprotion.tsv',sep = '\t')
  
  
  #Step2:Association analysis
  if (variables == 'continuous'){
    library(psych)
    
    cell_abundance = t(cell_abundance) %>% data.frame(check.names = F)
    
    phenotype_target = tibble::column_to_rownames(phenotype,sample_ID)
    phenotype_target =  phenotype_target[rownames(cell_abundance),]
    
    pearson <- corr.test(cell_abundance,phenotype_target,use="pairwise", method="pearson", 
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
    
    results = pearson_cor
    
    write.table(results,str_glue('iCWAS_result.tsv'),sep = '\t',row.names = F)
    
  } else if (variables == 'discrete'){
    
    cell_abundance_df = melt(cell_abundance)
    names(cell_abundance_df) = c('cell',sample_ID,'Abundance')
    merge_data <- left_join(cell_abundance_df,phenotype,by=sample_ID)
    
    groups = unique(merge_data[[sample_group]])
    
    if (length(groups) != 2) {
      stop("There must be exactly two groups in the sample group column for comparison.")
    }
    
    # Calculate log2FC and p-values
    results <- merge_data %>%
      group_by(cell) %>%
      summarise(
        mean_group1 = mean(Abundance[get(sample_group) == groups[1]], na.rm = TRUE),
        mean_group2 = mean(Abundance[get(sample_group) == groups[2]], na.rm = TRUE),
        log2FC = log2(mean_group2 + 1) - log2(mean_group1 + 1),
        p.value = t.test(Abundance ~ get(sample_group))$p.value
      )
    
    results <- results %>% 
      mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>% 
      mutate(p.adj = p.adjust(p.value, method = "BH"))
    
    g1 = groups[1]
    g2 = groups[2]
    names(results) = c("cell",paste0(g1,'_exp'),paste0(g2,'_exp'),
                       paste0("log2FC_",g2,'.',g1),"p.value","p.adj")
    
    write.table(results,str_glue('iCWAS_result.tsv'),sep = '\t',row.names = F)
    
  }
  
  return(results)
}
