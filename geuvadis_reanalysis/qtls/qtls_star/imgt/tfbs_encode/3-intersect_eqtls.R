devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- filter(gencode_chr_gene, gene_name %in% hla_genes) %>%
    select(gene_name, gene_id)

qtls <- 
    read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1) %>%
    select(gene_name, var_id, rank, var_from)

tfbs <- read_tsv("./chr6_tfbs_hg38.bed", col_names = FALSE)

join_list <- vector("list", nrow(qtls))
names(join_list) <- qtls$var_id

for(i in 1:nrow(qtls))
    join_list[[i]] <- 
	filter(tfbs, qtls$var_from[i] >= X2 & qtls$var_from[i] <= X3)

tfbs_qtls_intersect <- bind_rows(join_list, .id = "var_id") %>%
    select(var_id, tf = X4) %>%
    group_by(var_id) %>%
    summarize(tfbs = paste(tf, collapse = "/")) %>%
    ungroup()

left_join(qtls, tfbs_qtls_intersect, by = "var_id") %>%
    arrange(gene_name, rank) %>%
    write_tsv("./results.tsv")

