devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

whole_catalog <-
  read_tsv("./gwas_catalog_v1.0.1-associations_e86_r2016-11-28.tsv", guess_max = 20000) %>%
  select(LINK, `REPORTED GENE(S)`, `MAPPED_GENE`, `SNP_ID_CURRENT`, `MAPPED_TRAIT`) %>%
  mutate(SNP_ID_CURRENT = paste0("rs", SNP_ID_CURRENT))

catalog <- 
  read_tsv("./gwas_catalog_filtered.txt", col_names = c("variant", "trait"))

qtls <-
  read_qtltools("../../3-conditional_analysis/conditional_60_all.txt.gz") %>%
  filter(bwd_best == 1) %>%
  inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
  select(gene = gene_name, variant = var_id, rank)

rtc <- 
  list.files(".", pattern = "^rtc_results") %>%
  map_df(read_qtltools_rtc) %>%
  filter(rtc > .95) %>%
  inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
  select(gene = gene_name, qtl_var, gwas_var, d_prime, rtc)

qtls_rtc <-
  left_join(qtls, rtc, by = c("gene", "variant" = "qtl_var")) %>%
  left_join(catalog, by = c("gwas_var" = "variant")) %>%
  arrange(gene, rank, desc(rtc)) %>%
  mutate(rtc = round(rtc, 3)) %>%
  separate_rows(trait, sep = ";") %>%
  left_join(whole_catalog, by = c("gwas_var" = "SNP_ID_CURRENT")) %>% 
  group_by(gene, variant, rank, gwas_var, d_prime, rtc, trait) %>%
  summarise(link = paste(unique(LINK), collapse = " ")) %>%
  ungroup() %>%
  arrange(gene, rank)

write_tsv(qtls_rtc, "./results.tsv")
