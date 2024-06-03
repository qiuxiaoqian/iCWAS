setwd('/gold160t/qiuxiaoqian/SARS_PBMC/Figures/github_upload')

library(tidyr)
library(data.table)
library(psych)
library(ggpubr)



ELISPLOT_res = openxlsx::read.xlsx('data/ELISPOT.xlsx',rowNames = T,check.names = F)
log_ELISPLOT_res = log2(ELISPLOT_res+1)


cell_abundance=read.table('data/bulk_estimated_cell_proprotion.tsv',sep = '\t')
cell_abundance = t(cell_abundance) %>% data.frame(check.names = F)


for (virus in names(log_ELISPLOT_res)){
  
  print(virus)
  tmp = dplyr::select(log_ELISPLOT_res,virus)
  da = merge(tmp,cell_abundance,by='row.names')
  da = da[,-1]
  da = melt(da,id.vars = virus,variable.name = "cell_cluster",value.name = 'cell_proportion')
  da$cell_cluster = factor(da$cell_cluster,levels = sort(as.character(unique(da$cell_cluster))))
  
  p=ggplot(da, aes(x = cell_proportion,y = !!as.name(virus))) +
    theme_bw() +
    facet_wrap( . ~ cell_cluster,ncol = 7,scales = 'free')+
    labs(title = virus, x = 'Cell abundance (%)', y = bquote(log[2] ~ SFU/10^6 ~ cells)) +
    geom_point(color = 'black') +
    geom_smooth(method = 'lm', color = 'red', se = FALSE) +
    theme(legend.position = 'none') +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'pearson', 
             label.x.npc = 'left', label.y.npc = 'top', size = 2.7)+
    theme(
      plot.title = element_text(size = 12,color = "black",face = 'bold'),
      axis.text = element_text(size=10,color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10,colour = "black"),#angle = 60,hjust = 0.5, vjust = 0.3
      strip.background = element_rect(fill = NA,color=NA),
      legend.title.align=0.5)
  
  if (virus == "BA.5/BF.7"){virus = "BA.5.BF.7"}
  ggsave(str_glue('plot/iCWAS_correlation/{virus}.pdf'),p,width = 15,height = 20)
  
}
