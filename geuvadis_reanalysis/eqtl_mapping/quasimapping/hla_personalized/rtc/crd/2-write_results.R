devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

qtls <- 
    read_qtltools("../../3-conditional_analysis/conditional_50_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, variant = var_id, rank)

catalog <- read_tsv("./variants.tsv", col_names = c("var_id", "info"))

rtc_files <- list.files(".", pattern = "^rtc_results")  

rtc_crd <- rtc_files %>%
    map_df(read_qtltools_rtc) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, crd_var = gwas_var, d_prime, rtc)

qtls_rtc <- 
    left_join(qtls, rtc_crd, by = c("gene", "variant" = "qtl_var")) %>%
    left_join(catalog, by = c("crd_var" =  "var_id"), 
	      suffix = c("_imgt", "_crd")) %>%
    group_by(gene, variant, rank) %>%
    filter(rtc == max(rtc)) %>%
    ungroup()

write_tsv(qtls_rtc, "./results.tsv")
unlink(rtc_files)
