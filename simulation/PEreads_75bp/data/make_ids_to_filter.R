devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

hla_tags <- paste0("IMGT_", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

pri_ids <- 
    gencode_pri_tx %>%
    filter(gene_name %in% hla_genes) %>%
    pull(tx_id)

write(c(hla_tags, pri_ids), "./ids_to_filter.txt")

