devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_name, gene_id)

qtls <-
    read_tsv("../../2-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(gene, variant = var_id, rank)

catalog <-
    "~/hlaexpression/geuvadis_reanalysis/data/previous_qtls/top_qtl_catalog.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    separate_rows(X2, sep = ";") %>%
    separate(X2, c("study", "gene", "pvalue"), sep = "\\|") %>%
    mutate(pvalue = round(as.numeric(pvalue), 2)) %>%
    select(variant = X1, gene, study, pvalue)

rtc_files <- "./rtc_results.txt" 

rtc <- read_qtltools_rtc(rtc_files) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_ref = gwas_var, d_prime, rtc)

qtls_rtc <-
    left_join(qtls, rtc, by = c("gene", "variant" =  "qtl_var")) %>%
    left_join(catalog, by = c("qtl_ref" =  "variant"), 
	      suffix = c("_personalized", "_ref")) %>%
    filter(gene_personalized == gene_ref) %>%
    group_by(gene_personalized, rank) %>%
    filter(rtc == max(rtc)) %>%
    ungroup() %>%
    mutate(rtc = round(rtc, 2),
	   info = paste0(study, " (", pvalue, ")")) %>%
    select(gene = gene_personalized, rank, eqtl_personalized = variant, 
	   eqtl_previous = qtl_ref, d_prime, rtc, info) %>%
    arrange(gene, rank)

write_tsv(qtls_rtc, "./results.tsv")
unlink(rtc_files)
