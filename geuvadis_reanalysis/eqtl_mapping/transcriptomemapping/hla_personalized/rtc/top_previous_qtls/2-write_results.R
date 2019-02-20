devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_name, gene_id)

qtls <- read_tsv("../../2-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(gene, rank, qtl_personalized = var_id)

catalog <-
    "~/hlaexpression/geuvadis_reanalysis/data/previous_qtls/top_qtl_catalog.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    separate_rows(X2, sep = ";") %>%
    separate(X2, c("study", "gene", "pvalue"), sep = "\\|") %>%
    mutate(pvalue = round(as.numeric(pvalue), 2)) %>%
    select(variant = X1, gene, study, pvalue)

rtc <- read_qtltools_rtc("./rtc_results.txt") %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_ref = gwas_var, r2 = r_squared, d_prime, rtc)

qtls_rtc <-
    left_join(qtls, rtc, by = c("gene", "qtl_personalized" =  "qtl_var")) %>%
    left_join(catalog, by = c("qtl_ref" =  "variant"), 
	      suffix = c("_personalized", "_ref")) %>%
    filter(gene_personalized == gene_ref) %>% 
    mutate(rtc = round(rtc, 2),
	   r2 = round(r2, 2),
	   d_prime = round(d_prime, 2),
	   info = paste0(study, " (", pvalue, ")")) %>%
    select(gene = gene_personalized, rank, qtl_personalized, qtl_previous = qtl_ref, 
	   r2, d_prime, rtc, info)

write_tsv(qtls_rtc, "./results.tsv")
