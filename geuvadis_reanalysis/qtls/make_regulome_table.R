devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

imgt_qtl <- 
    "./qtls_star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(gene_name, rank, var_id)

pri_qtl <-
    "./qtls_star/pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(gene_name, rank, var_id)

qtl_df <-
    list(imgt = imgt_qtl, pri = pri_qtl) %>%
    bind_rows(.id = "index") %>%
    select(gene_name, rank, index, var_id) %>%
    arrange(gene_name, rank, index)

regulome <- read_tsv("./regulomedb_results.bed", col_names = FALSE) %>%
    select(X4) %>%
    separate(X4, c("var_id", "score"), sep = ";") %>%
    distinct()

qtl_df <- left_join(qtl_df, regulome)

write_tsv(qtl_df, "./qtls_with_regulome_score.tsv")
