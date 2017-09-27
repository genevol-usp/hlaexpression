devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

pri_ids <- 
    gencode_pri_tx %>%
    filter(gene_name %in% hla_genes) %>%
    pull(tx_id)

imgt_ids <- 
    read_tsv("./genos.tsv") %>%
    pull(allele) %>%
    unique() %>%
    sort()

write(c(imgt_ids, pri_ids), "./ids_to_filter.txt")

