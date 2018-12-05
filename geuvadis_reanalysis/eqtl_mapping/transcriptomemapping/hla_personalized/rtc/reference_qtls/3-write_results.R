devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_name, gene_id)

qtls <-
    read_tsv("../../2-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(gene, variant = var_id, rank)

ref_qtls <-
    read_tsv("../../../reference/3-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(gene, variant = var_id, rank)

rtc_files <- "./rtc_results.txt" 

rtc <- read_qtltools_rtc(rtc_files) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_ref = gwas_var, d_prime, rtc)

qtls_rtc <-
    left_join(qtls, rtc, by = c("gene", "variant" =  "qtl_var")) %>%
    left_join(ref_qtls, by = c("qtl_ref" =  "variant"), 
	      suffix = c("_personalized", "_ref")) %>%
    filter(gene_personalized == gene_ref) %>%
    group_by(gene_personalized, rank_personalized) %>%
    filter(rtc == max(rtc)) %>%
    ungroup() %>%
    select(gene = gene_personalized, variant_personalized = variant, 
	   rank_personalized, variant_ref = qtl_ref, rank_ref, d_prime, rtc) %>%
    mutate(rtc = round(rtc, 2))

write_tsv(qtls_rtc, "./results.tsv")
unlink(rtc_files)
