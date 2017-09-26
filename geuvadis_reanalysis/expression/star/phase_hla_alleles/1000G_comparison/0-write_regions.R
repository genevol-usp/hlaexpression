devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

bed <-
  gencode_chr_gene %>%
  filter(gene_name %in% hla_genes) %>%
  select(chr, start, end)

write_tsv(bed, "./hla_regions.bed", col_names = FALSE)
