devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- filter(gencode_chr_gene, gene_name %in% hla_genes)

eqtl <- 
    "../../../qtls/qtls_star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, rank == 0L, bwd_best == 1L) %>%
    select(phen_id, var_id, var_chr, var_from, var_to) 

write_lines(eqtl$var_id, "./best_eqtl_rsID.txt")
write_tsv(eqtl, "./best_eqtl.tsv")
