devtools::load_all("~/hlaseqlib")
library(tidyverse)

imgt_qtls <-
    read_qtltools("../3-conditional_analysis/conditional_70_all.txt.gz") %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(gene_id = phen_id, gene_name, rank, var_id, pos = var_from)

pri_qtls <-
    read_qtltools("../../pri/3-conditional_analysis/conditional_60_all.txt.gz") %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(gene_id = phen_id, gene_name, rank, var_id, pos = var_from)

results_best_imgt <- read_tsv("./results.best") %>%
    separate(GENE, c("gene_id", "chr", "pos", "ref", "alt"), sep = "_", convert = TRUE) %>%
    inner_join(imgt_qtls, by = c("gene_id", "pos")) %>%
    select(gene_name, rank, var_id, pos, caveman_var_pos = POS, P, CaVEMaN, Probability) %>%
    arrange(gene_name, rank) %>%
    mutate(P = -log10(P),
	   best = as.integer(pos == caveman_var_pos)) %>%
    group_by(gene_name, rank) %>%
    filter((n() > 1 & any(best == 1) & best == 1) | (all(best == 0) | all(best == 1))) %>% 
    group_by(gene_name, rank, var_id, pos, P, CaVEMaN, Probability, best) %>%
    summarize(caveman_var_pos = paste(caveman_var_pos, collapse = "/")) %>%
    ungroup()

results_best_pri <- read_tsv("../../pri/caveman/results.best") %>%
    separate(GENE, c("gene_id", "chr", "pos", "ref", "alt"), sep = "_", convert = TRUE) %>%
    inner_join(pri_qtls, by = c("gene_id", "pos")) %>%
    select(gene_name, rank, var_id, pos, caveman_var_pos = POS, P, CaVEMaN, Probability) %>%
    arrange(gene_name, rank) %>%
    mutate(P = -log10(P),
	   best = as.integer(pos == caveman_var_pos)) %>%
    group_by(gene_name, rank) %>%
    filter((n() > 1 & any(best == 1) & best == 1) | (all(best == 0) | all(best == 1))) %>% 
    group_by(gene_name, rank, var_id, pos, P, CaVEMaN, Probability, best) %>%
    summarize(caveman_var_pos = paste(caveman_var_pos, collapse = "/")) %>%
    ungroup()

print(results_best_imgt, n = Inf)
print(results_best_pri, n = Inf)

list(HLA_personalized = results_best_imgt, Reference = results_best_pri) %>%
    bind_rows(.id = "index") %>%
    write_tsv("./results.hla")

