devtools::load_all("~/hlaseqlib")
library(tidyverse)

pri_ids <- 
    gencode_pri_tx %>%
    filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
    pull(tx_id)

imgt_ids <- 
    read_tsv("../../data/genos.tsv") %>%
    pull(allele) %>%
    unique()

write(c(imgt_ids, pri_ids), "./ids_to_filter.txt")

