devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

pri_ids <- 
    gencode_pri_tx %>%
    filter(gene_name %in% hla_genes) %>%
    pull(tx_id)

write(pri_ids, "./hla_ids_pri.txt")
