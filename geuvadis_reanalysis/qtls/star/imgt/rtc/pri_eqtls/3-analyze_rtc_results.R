devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

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

rtc_pri <- 
    read_qtltools_rtc("./rtc_results.txt") %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_pri = gwas_var, d_prime, rtc)

qtls_rtc_pri <- 
    left_join(qtls, rtc_pri, by = c("gene", "variant" =  "qtl_var")) %>%
    left_join(pri_eqtls, by = c("qtl_pri" =  "variant"), 
	      suffix = c("_imgt", "_pri")) %>%
    drop_na() %>%
    filter(rtc > 0.95) %>%
    select(gene_imgt, variant_imgt = variant, rank_imgt,
	   gene_pri, variant_pri = qtl_pri, rank_pri, d_prime, rtc) %>%
    mutate(rtc = round(rtc, 3)) %>%
    arrange(gene_imgt, rank_imgt, rank_pri)

write_tsv(qtls_rtc_pri, "./results.tsv")
