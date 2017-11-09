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

previous_eqtls <-
    read_tsv("./previous_eQTL.tsv", col_names = c("variant", "study"))

rtc_previous <- 
    read_qtltools_rtc("./rtc_results.txt") %>% 
    rename(qtl_previous = gwas_var) %>%
    inner_join(gencode_hla, by = c("gene" = "gene_id")) %>%
    select(gene = gene_name, qtl_var, qtl_previous, d_prime, rtc)

qtls_rtc_previous <- 
    left_join(qtls, rtc_previous, 
	      by = c("gene", "variant" =  "qtl_var")) %>%
    inner_join(previous_eqtls, by = c("qtl_previous" =  "variant")) %>%
    filter(rtc > 0.9) %>%
    mutate(rtc = round(rtc, 3)) %>%
    distinct() %>%
    group_by(gene, variant, rank, qtl_previous, d_prime, rtc) %>%
    summarise(study = paste(study, collapse = "/")) %>%
    ungroup() %>%
    arrange(gene, rank, desc(rtc))

write_tsv(qtls_rtc_previous, "./results.tsv")
