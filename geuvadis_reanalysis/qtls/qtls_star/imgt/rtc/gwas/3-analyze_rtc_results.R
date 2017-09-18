devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_id, gene_name)

whole_catalog <-
  read_tsv("./gwas_catalog_v1.0.1-associations_e86_r2016-11-28.tsv", guess_max = 20000) %>%
  select(LINK, `REPORTED GENE(S)`, `MAPPED_GENE`, `SNP_ID_CURRENT`, `MAPPED_TRAIT`) %>%
  mutate(SNP_ID_CURRENT = paste0("rs", SNP_ID_CURRENT))

catalog <- 
  read_tsv("./gwas_catalog_filtered.txt", col_names = c("variant", "trait"))

qtls <-
  read_qtltools("../../conditional_analysis/conditional_60_all.txt.gz") %>%
  filter(bwd_best == 1) %>%
  inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
  select(gene = gene_name, variant = var_id, rank)

rtc_names <- 
  c("gwas_variant", "qtl_variant", "gene", "gene_group", "gwas_variant_chr",
    "gwas_variant_pos", "gwas_variant_rank", "qtl_variant_chr", 
    "qtl_variant_pos", "qtl_variant_rank", "gene_chr", "gene_pos", 
    "distance_variants", "distance_gwas_gene", "gwas_variant_region_index",
    "qtl_variant_region_index", "region_start", "region_end", "n_variants",
    "rtc", "d_prime", "r_squared")

rtc <- 
  read_delim("./rtc_results.txt", delim = " ", col_names = rtc_names) %>% 
  filter(rtc > .9) %>%
  inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
  select(gene = gene_name, qtl_variant, gwas_variant, distance_variants, distance_gwas_gene, rtc)

qtls_rtc <-
  left_join(qtls, rtc, by = c("gene", "variant" = "qtl_variant")) %>%
  group_by(gene, rank) %>%
  filter(rtc == max(rtc)) %>%
  ungroup() %>%
  left_join(catalog, by = c("gwas_variant" = "variant")) %>%
  arrange(gene, rank, desc(rtc), distance_gwas_gene) %>%
  mutate(rtc = round(rtc, 3)) %>%
  left_join(whole_catalog, by = c("gwas_variant" = "SNP_ID_CURRENT")) %>% 
  group_by(gene, variant, rank, gwas_variant, distance_variants, 
           distance_gwas_gene, rtc, trait) %>%
  summarise(link = paste(unique(LINK), collapse = " ")) %>%
  ungroup() %>%
  arrange(gene, rank)

write_tsv(qtls_rtc, "./results.tsv")
