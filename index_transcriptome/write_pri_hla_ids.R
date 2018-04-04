devtools::load_all("~/hlaseqlib")
library(tidyverse)

imgt_loci <- readLines("../imgt_index_v2/imgt_loci.txt") %>% 
    paste0("HLA-", .)

pri_ids <- gencode_pri_tx %>%
    filter(gene_name %in% imgt_loci) %>%
    pull(tx_id)

write(pri_ids, "./hla_ids_pri.txt")

