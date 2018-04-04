library(tidyverse)

qtl_table <- 
    read_tsv("../../../reference/3-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1L) %>%
    mutate(info = paste(gene, "PRI", sep = ": ")) %>%
    select(var_id, info)

write_tsv(qtl_table, "./ref_qtls.tsv", col_names = FALSE)
