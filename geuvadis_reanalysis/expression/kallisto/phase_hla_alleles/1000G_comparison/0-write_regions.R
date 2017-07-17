devtools::load_all("~/genomicRutils")
library(tidyverse)

bed <-
  gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(chr, start, end)

write_tsv(bed, "./hla_regions.bed", col_names = FALSE)
