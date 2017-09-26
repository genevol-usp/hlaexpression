devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))

gencode_hla <- filter(gencode_chr_gene, gene_name %in% hla_genes)

qtls <-
  read_qtltools("../../conditional_analysis/conditional_60_all.txt.gz") %>%
  filter(bwd_best == 1) %>%
  inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
  select(gene = gene_name, variant = var_id, rank)

previous_eqtls <- read_tsv("./catalog.tsv", col_names = c("variant", "info"))

rtc_names <- 
  c("qtl_previous", "qtl_variant", "gene", "gene_group", "gwas_variant_chr",
    "gwas_variant_pos", "gwas_variant_rank", "qtl_variant_chr", 
    "qtl_variant_pos", "qtl_variant_rank", "gene_chr", "gene_pos", 
    "distance_variants", "distance_eqtl_previous_gene", "gwas_variant_region_index",
    "qtl_variant_region_index", "region_start", "region_end", "n_variants",
    "rtc", "d_prime", "r_squared")

rtc_previous <- 
  read_delim("./rtc_GeuvadisBestQTLs.txt", delim = " ", 
             col_names = rtc_names) %>% 
  inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
  select(gene = gene_name, qtl_variant, qtl_previous, distance_variants,
         distance_eqtl_previous_gene, rtc)

qtls_rtc_previous <- 
  left_join(qtls, rtc_previous, by = c("gene", "variant" =  "qtl_variant")) %>%
  inner_join(previous_eqtls, by = c("qtl_previous" =  "variant")) %>%
  group_by(gene, rank) %>%
  filter(rtc == max(rtc)) %>%
  ungroup() %>%
  mutate(rtc = round(rtc, 3)) %>%
  distinct() %>%
  separate(info, c("geuvadis_variant", "geuvadis_gene", "geuvadis_eff_size", "geuvadis_pvalue"), sep = ":") %>% 
  group_by(gene, variant, rank, geuvadis_variant, rtc) %>%
  summarize(geuvadis_gene = paste(geuvadis_gene, collapse = "/")) %>%
  select(gene, variant, rank, geuvadis_gene, geuvadis_variant, rtc) %>% 
  arrange(gene, rank)

write_tsv(qtls_rtc_previous, "./results.tsv")
