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

catalog <- read_tsv("./variants.tsv", col_names = c("var_id", "info"))

rtc <- 
    read_qtltools_rtc("./rtc_results.txt") %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, rtc_var = gwas_var, d_prime, rtc)

qtls_rtc <- 
    left_join(qtls, rtc, by = c("gene", "variant" =  "qtl_var")) %>%
    left_join(catalog, by = c("rtc_var" =  "var_id"), 
	      suffix = c("_imgt", "_rtc")) %>%
    drop_na() %>%
    filter(rtc > 0.9) %>%
    mutate_at(vars(rtc, d_prime), ~round(., 3)) %>%
    arrange(gene, rank, desc(rtc))

write_tsv(qtls_rtc, "./results.tsv")
