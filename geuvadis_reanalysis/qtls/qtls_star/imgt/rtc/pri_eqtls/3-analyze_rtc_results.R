devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

qtls <-
    read_qtltools("../../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, variant = var_id, rank)

pri_eqtls <- 
    "../../../pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
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

rtc_pri <- 
    read_delim("./rtc_results.txt", delim = " ", col_names = rtc_names) %>% 
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_variant, qtl_pri = gwas_variant, rtc)

qtls_rtc_pri <- 
    left_join(qtls, rtc_pri, by = c("gene", "variant" =  "qtl_variant")) %>%
    left_join(pri_eqtls, by = c("qtl_pri" =  "variant"), 
	      suffix = c("_imgt", "_pri")) %>%
    drop_na() %>%
    filter(rtc > 0.9) %>%
    select(gene_imgt, variant_imgt = variant, rank_imgt,
	   gene_pri, variant_pri = qtl_pri, rank_pri, rtc) %>%
    mutate(rtc = round(rtc, 3)) %>%
    arrange(gene_imgt, rank_imgt, rank_pri)

write_tsv(qtls_rtc_pri, "./results.tsv")
