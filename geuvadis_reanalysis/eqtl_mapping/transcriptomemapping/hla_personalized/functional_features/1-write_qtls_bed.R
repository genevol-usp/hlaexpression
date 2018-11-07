library(tidyverse)

qtls <- read_tsv("../2-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(var_chr, var_from, var_to, gene, var_id, rank) %>%
    mutate(var_chr = paste0("chr", var_chr),
	   var_from = var_from - 1L) %>%
    unite(info, gene:rank, sep = ":")

write_tsv(qtls, "./hla.qtls.bed", col_names = FALSE)
