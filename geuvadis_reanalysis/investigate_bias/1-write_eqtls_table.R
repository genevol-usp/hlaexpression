devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

eqtl_imgt <- 
    "../eqtl_mapping/transcriptomemapping/hla_personalized/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1L) %>%
    group_by(phen_id) %>%
    slice(which.min(bwd_pval)) %>%
    ungroup() %>%
    select(phen_id, var_id, bwd_slope) 

eqtl_conventional <- 
    "../eqtl_mapping/genomemapping/3-conditional_analysis/conditional_70_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1L) %>%
    group_by(phen_id) %>%
    slice(which.min(bwd_pval)) %>%
    ungroup() %>%
    select(phen_id, var_id, bwd_slope) 

eqtl_df <- list(personalized = eqtl_imgt, conventional = eqtl_conventional) %>%
    bind_rows(.id = "index")

write_tsv(eqtl_df, "./best_eqtl.tsv")
write_lines(unique(eqtl_df$var_id), "./best_eqtl_rsID.txt")
