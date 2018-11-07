devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

catalog <- 
    "~/hlaexpression/geuvadis_reanalysis/data/gwas_catalog/gwas_catalog_filtered.txt" %>% 
    read_tsv(col_names = c("variant", "trait"))

qtls <-
    read_qtltools("../../2-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, variant = var_id, rank)

rtc_files <- list.files(".", pattern = "^rtc_results")   

rtc <- rtc_files %>%
    map_df(read_qtltools_rtc) %>%
    filter(rtc > .95) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, gwas_var, d_prime, rtc)

qtls_rtc <- left_join(qtls, rtc, by = c("gene", "variant" = "qtl_var")) %>%
    left_join(catalog, by = c("gwas_var" = "variant")) %>%
    separate_rows(trait, sep = ";") %>%
    group_by(gene, variant, rank, trait) %>%
    summarise(gwas_var = paste(gwas_var, collapse = ", ")) %>%
    ungroup() %>%
    mutate(`trait (GWAS variant)` = paste0(trait, " (", gwas_var, ")")) %>%
    select(gene, rank, variant, `trait (GWAS variant)`) %>%
    group_by(gene, rank, variant) %>%
    summarise(`trait (GWAS variant)` = paste(`trait (GWAS variant)`, collapse = "; ")) %>%
    ungroup() %>%
    mutate(gene = factor(gene, levels = gencode_hla$gene_name)) %>%
    arrange(gene, rank)

write_tsv(qtls_rtc, "./results.tsv")
unlink(rtc_files)
