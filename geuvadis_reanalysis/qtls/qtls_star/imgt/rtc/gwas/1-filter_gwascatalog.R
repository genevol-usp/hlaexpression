library(tidyverse)

vcf <- 
    "zcat ../../../../genotypes/eur_maf05.vcf.gz | grep -v '##' | cut -f1,2,3" %>%
    data.table::fread() %>%
    as.tibble() 

gwascat <-
    read_tsv("./gwas_catalog_v1.0.1-associations_e86_r2016-11-28.tsv", 
	     guess_max = 1e5) %>% 
    filter(`P-VALUE` <= 1e-8) %>%
    select(ID = SNP_ID_CURRENT, trait = MAPPED_TRAIT) %>%
    mutate(ID = paste0("rs", ID)) %>%
    inner_join(vcf, by = "ID") %>%
    arrange(`#CHROM`, POS) %>%
    select(ID, trait)

write_tsv(gwascat, "./gwas_catalog_filtered.txt", col_names = FALSE)
