devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

eqtl_imgt <- 
    "../qtls/star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1L) %>%
    group_by(phen_id) %>%
    slice(which.min(bwd_pval)) %>%
    ungroup() %>%
    select(phen_id, var_id, var_chr, var_from, var_to) 

eqtl_pri <- 
    "../qtls/star/pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1L) %>%
    group_by(phen_id) %>%
    slice(which.min(bwd_pval)) %>%
    ungroup() %>%
    select(phen_id, var_id, var_chr, var_from, var_to) 

eqtl_df <- list(imgt = eqtl_imgt, ref = eqtl_pri) %>%
    bind_rows(.id = "index")

write_tsv(eqtl_df, "./best_eqtl.tsv")
write_lines(unique(eqtl_df$var_id), "./best_eqtl_rsID.txt")