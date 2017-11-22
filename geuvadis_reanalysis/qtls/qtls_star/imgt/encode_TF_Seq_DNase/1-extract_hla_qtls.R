devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_name, gene_id)

qtls <-
    read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1) %>%
    select(var_chr, var_from, var_to, gene_name, var_id, rank) %>%
    mutate(var_chr = paste0("chr", var_chr),
	   var_from = var_from - 1L) %>%
    unite(info, gene_name:rank, sep = ":")

write_tsv(qtls, "./hla.qtls.bed", col_names = FALSE)
