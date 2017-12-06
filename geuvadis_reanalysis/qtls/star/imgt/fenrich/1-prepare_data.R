devtools::load_all("~/hlaseqlib")
library(tidyverse)

qtls <- read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz")

bed_significant <- qtls %>% 
    filter(phen_id %in% gencode_hla$gene_id, rank == 0L, bwd_best == 1L, 
	   bwd_signif == 1L) %>%
    select(var_chr, var_from, var_to, var_id, phen_id, strand) %>%
    mutate(var_chr = as.integer(var_chr),
	   var_from = var_from - 1L) %>%
    arrange(var_chr, var_from)
  
write_tsv(bed_significant, "./genes_significant.bed", col_names = FALSE)

bed_quantified <-
    read_delim("../2-permutations/results/permutations_60.txt.gz", delim = " ", 
	       col_names = FALSE) %>%
    select(X2, X3, X4, X1, X8, X5) %>%
    mutate(X3 = X3 - 1L) %>%
    arrange(X2, X3)

write_tsv(bed_quantified, "./genes_quantified.bed", col_names = FALSE)

read_tsv("../tfbs_encode/tfbs_hg38.bed", col_names = FALSE) %>%
    filter(X1 %in% paste0("chr", 1:22)) %>%
    mutate(X1 = as.integer(sub("^chr", "", X1))) %>%
    arrange(X1, X2) %>%
    write_tsv("./tfbs_encode_hg38.bed", col_names = FALSE)
    
