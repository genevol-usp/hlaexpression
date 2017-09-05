devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_pri_tx %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  pull(tx_id) %>%
  c(paste0("IMGT_", c("A", "B", "C", "DQA1", "DQB1", "DRB1")), .) %>%
  writeLines("./ids_to_filter.txt")

