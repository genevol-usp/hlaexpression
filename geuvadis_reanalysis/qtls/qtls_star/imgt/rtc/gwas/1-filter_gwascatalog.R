library(tidyverse)

gwascat <-
    read_tsv("./gwas_catalog_v1.0.1-associations_e86_r2016-11-28.tsv", 
	     guess_max = 1e5) %>% 
    filter(!is.na(SNP_ID_CURRENT), `P-VALUE` <= 1e-8) %>%
    select(ID = SNP_ID_CURRENT, trait = MAPPED_TRAIT) %>%
    mutate(ID = paste0("rs", ID)) %>%
    select(ID, trait) %>%
    group_by(ID) %>%
    summarize(trait = paste(unique(trait), collapse = ";")) %>%
    ungroup()

write_tsv(gwascat, "./gwas_catalog_filtered.txt", col_names = FALSE)
