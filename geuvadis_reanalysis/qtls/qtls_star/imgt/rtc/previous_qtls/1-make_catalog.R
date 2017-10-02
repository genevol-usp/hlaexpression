library(tidyverse)

qtl_table <-
  read_tsv("../../../../../data/previous_qtls/previous_eQTL.tsv") %>%
  mutate(study = paste(paste0(phenotype, ":"), source)) %>%
  select(variant, study)

write_tsv(qtl_table, "./previous_eQTL.tsv", col_names = FALSE)
