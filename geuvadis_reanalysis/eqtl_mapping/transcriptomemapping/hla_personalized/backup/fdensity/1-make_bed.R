devtools::load_all("~/hlaseqlib")
library(tidyverse)

qtl <- read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1) %>%
    select(var_chr, var_from, var_to, var_id, phen_id, strand) %>%
    mutate(var_chr = paste0("chr", var_chr),
	   var_from = var_from - 1L)

write_tsv(qtl, "./hla_qtls.bed", col_names = FALSE)
