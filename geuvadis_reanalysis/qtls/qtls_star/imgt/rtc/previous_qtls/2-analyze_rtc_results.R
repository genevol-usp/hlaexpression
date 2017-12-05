devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

qtls <-
    read_qtltools("../../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rsid = var_id, rank)

previous_eqtls <-
    "/home/vitor/hlaexpression/geuvadis_reanalysis/data/previous_qtls/qtl_catalog.tsv" %>%
    read_tsv(col_names = c("rsid", "info")) %>%
    separate_rows(info, sep = ";") %>%
    separate(info, c("study", "tissue", "gene", "pvalue"), sep = "\\|")
    
rtc <- 
    list.files(".", pattern = "^rtc_results") %>%
    map_df(read_qtltools_rtc) %>%
    rename(qtl_previous = gwas_var) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_previous, d_prime, rtc)

rtc_df <- 
    left_join(qtls, rtc, by = c("gene", "rsid" =  "qtl_var")) %>%
    left_join(previous_eqtls, by = c("gene", "qtl_previous" =  "rsid")) %>%
    filter(rtc > 0.95, !is.na(study)) %>%
    group_by(gene, rsid) %>%
    filter(rtc == max(rtc)) %>%
    group_by(gene, rsid, rank, qtl_previous, d_prime, rtc, study) %>%
    summarize(tissue = paste(tissue, collapse = ",")) %>%
    ungroup() %>%
    arrange(gene, rank, desc(rtc)) %>%
    mutate(rtc = round(rtc, 3))

write_tsv(rtc_df, "./results.tsv")
