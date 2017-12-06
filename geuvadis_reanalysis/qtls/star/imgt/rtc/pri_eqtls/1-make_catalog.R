devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

hla_ids <- filter(gencode_chr_gene, gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

qtl_table <-
    "../../../pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_best == 1L) %>%
    inner_join(hla_ids, by = c("phen_id" = "gene_id")) %>%
    mutate(info = paste(gene_name, "PRI", sep = ": ")) %>%
    select(var_id, info)

write_tsv(qtl_table, "./pri_eQTL.tsv", col_names = FALSE)
