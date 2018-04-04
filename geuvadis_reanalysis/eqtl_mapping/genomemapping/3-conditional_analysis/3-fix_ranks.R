devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

qtltools_res <- 
    read_qtltools("./conditional_70_all.txt.gz") %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    mutate(tss = ifelse(strand == "+", phen_from, phen_to),
	   dist_to_tss = var_from - tss,
	   pval = -log10(bwd_pval)) %>%
    select(gene = gene_name, tss, strand, rank, var_chr, var_id, var_from, 
	   var_to, dist, dist_to_tss, pval, best = bwd_best, 
	   signif = bwd_signif) 

fix_rank <- qtltools_res %>%
    filter(best == 1) %>%
    arrange(gene, desc(pval), rank) %>%
    group_by(gene) %>%
    mutate(new_rank = 1:n() - 1L) %>%
    ungroup() %>%
    select(gene, rank, new_rank)

out_df <- left_join(qtltools_res, fix_rank, by = c("gene", "rank")) %>%
    select(gene, tss, strand, rank = new_rank, var_chr, var_id, var_from, var_to,
	   dist, dist_to_tss, pval, best, signif) %>%
    group_by(gene, var_id) %>%
    filter(pval == max(pval)) %>%
    ungroup() %>%
    mutate(gene = reorder(gene, tss), rank = as.character(rank))

write_tsv(out_df, "hla_qtls.tsv")
