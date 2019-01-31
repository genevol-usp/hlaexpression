devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_name, gene_id)

qtls <- read_tsv("../../../2-conditional_analysis/hla_qtls_all.tsv") 

best_qtls <- qtls %>%
    filter(best == 1) %>%
    select(gene, rank, qtl_hlapers = var_id)

catalog <-
    "/home/vitor/hlaexpression/geuvadis_reanalysis/data/previous_qtls/top_qtl_catalog.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    separate_rows(X2, sep = ";") %>%
    separate(X2, c("study", "gene", "pvalue"), sep = "\\|") %>%
    mutate(pvalue = round(as.numeric(pvalue), 2)) %>%
    select(variant = X1, gene, study, pvalue) %>%
    filter(study %in% c("petersdorf2015", "kulkarni2011", "vince2017", "xl9_raj2016", "ou2018"))

ref_vars_status <- catalog %>%
    select(-pvalue) %>%
    inner_join(qtls, by = c("gene", "variant" = "var_id")) %>%
    filter(signif == 1) %>%
    select(validated_var = variant, gene, study, rank, 
	   validated_var_pval = pval, validated_var_signif = signif)

rtc <- read_qtltools_rtc("../rtc_results.txt") %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_ref = gwas_var, r_squared, d_prime, rtc)

qtls_rtc <-
    left_join(best_qtls, rtc, by = c("gene", "qtl_hlapers" = "qtl_var")) %>%
    inner_join(catalog, by = c("qtl_ref" =  "variant"), 
	      suffix = c("_hlapers", "_ref")) %>% 
    filter(gene_hlapers == gene_ref)

af <- 
    "/home/vitor/hlaexpression/geuvadis_reanalysis/data/genotypes/validated_af.tsv" %>%
    read_delim(col_names = c("varid", "af"), delim = " ")

out <- qtls_rtc %>%
    select(gene = gene_hlapers, rank, qtl = qtl_hlapers, 
	   validated_var = qtl_ref, study, r_squared, d_prime, rtc) %>%
    left_join(ref_vars_status, by = c("gene", "rank", "validated_var", "study")) %>%
    left_join(af, by = c("qtl" = "varid")) %>% 
    left_join(af, by = c("validated_var" = "varid"), suffix = c("_qtl", "_validated")) %>% 
    mutate_at(vars(r_squared:validated_var_pval), ~round(., 2)) %>%
    select(gene, rank, qtl, validated_var, study, r_squared, d_prime, af_qtl, af_validated, everything())


write_csv(out, "./validated_variants_results.csv")
