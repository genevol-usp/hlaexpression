devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

hla_tags <- gsub("HLA-", "IMGT_", hla_genes)
    
pri_ids <- gencode_pri_tx %>%
    filter(gene_name %in% hla_genes) %>%
    pull(tx_id)

write(c(hla_tags, pri_ids), "./ids_to_filter.txt")
write(gencode_hla$gene_id, "./geneids_to_filter.txt")
