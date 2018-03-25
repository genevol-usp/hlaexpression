devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

qtl_table <-
    "../../../reference/3-conditional_analysis/conditional_50_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_best == 1L) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    mutate(info = paste(gene_name, "PRI", sep = ": ")) %>%
    select(var_id, info)

write_tsv(qtl_table, "./reference_qtls.tsv", col_names = FALSE)
