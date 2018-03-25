devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

qtls <-
    read_qtltools("../../3-conditional_analysis/conditional_50_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rsid = var_id, rank)

previous_eqtls <-
    "/home/vitor/hlaexpression/geuvadis_reanalysis/data/previous_qtls/qtl_catalog.tsv" %>%
    read_tsv(col_names = c("rsid", "info")) %>%
    separate_rows(info, sep = ";") %>%
    separate(info, c("study", "gene", "pvalue"), sep = "\\|")
   
rtc_files <- list.files(".", pattern = "^rtc_results")

rtc <- rtc_files %>%
    map_df(read_qtltools_rtc) %>%
    rename(qtl_previous = gwas_var) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_previous, d_prime, rtc)

rtc_df <- 
    left_join(qtls, rtc, by = c("gene", "rsid" =  "qtl_var")) %>%
    left_join(previous_eqtls, by = c("gene", "qtl_previous" =  "rsid")) %>%
    group_by(gene, rsid) %>%
    filter(rtc == max(rtc)) %>%
    ungroup() %>%
    mutate(rtc = round(rtc, 3),
	   pvalue = round(as.numeric(pvalue), 1),
	   study_pval = paste0(study, " (", pvalue, ")")) %>% 
    group_by(gene, rsid, rank, qtl_previous, d_prime, rtc) %>%
    summarise(study_pval = paste(study_pval, collapse = "/")) %>%
    ungroup() %>%
    left_join(gencode_hla, by = c("gene" = "gene_name")) %>%
    arrange(start, rank) %>%
    select(gene, rank, rsid, qtl_previous, d_prime, rtc, study_pval)
    
write_tsv(rtc_df, "./results.tsv")
unlink(rtc_files)
