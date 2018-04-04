devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

imgt_qtls <- read_tsv("../3-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1L) %>%
    left_join(gencode_hla, by = c("gene" = "gene_name")) %>%
    select(gene_id, gene, rank, var_id, pos = var_from)

ref_qtls <- read_tsv("../../reference/3-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1L) %>%
    left_join(gencode_hla, by = c("gene" = "gene_name")) %>%
    select(gene_id, gene, rank, var_id, pos = var_from)

results_best_imgt <- read_tsv("./results.best") %>%
    separate(GENE, c("gene_id", "chr", "pos", "ref", "alt"), sep = "_", convert = TRUE) %>%
    inner_join(imgt_qtls, by = c("gene_id", "pos")) %>%
    select(gene, rank, var_id, pos, caveman_var_pos = POS, P, CaVEMaN, Probability) %>%
    arrange(gene, rank) %>%
    mutate(P = -log10(P),
	   best = as.integer(pos == caveman_var_pos)) %>%
    group_by(gene, rank) %>%
    filter((n() > 1 & any(best == 1) & best == 1) | (all(best == 0) | all(best == 1))) %>% 
    group_by(gene, rank, var_id, pos, P, CaVEMaN, Probability, best) %>%
    summarize(caveman_var_pos = paste(caveman_var_pos, collapse = "/")) %>%
    ungroup()

results_best_ref <- read_tsv("../../reference/caveman/results.best") %>%
    separate(GENE, c("gene_id", "chr", "pos", "ref", "alt"), sep = "_", convert = TRUE) %>%
    inner_join(ref_qtls, by = c("gene_id", "pos")) %>%
    select(gene, rank, var_id, pos, caveman_var_pos = POS, P, CaVEMaN, Probability) %>%
    arrange(gene, rank) %>%
    mutate(P = -log10(P),
	   best = as.integer(pos == caveman_var_pos)) %>%
    group_by(gene, rank) %>%
    filter((n() > 1 & any(best == 1) & best == 1) | (all(best == 0) | all(best == 1))) %>% 
    group_by(gene, rank, var_id, pos, P, CaVEMaN, Probability, best) %>%
    summarize(caveman_var_pos = paste(caveman_var_pos, collapse = "/")) %>%
    ungroup()

list(HLA_personalized = results_best_imgt, Ref_Transcriptome = results_best_ref) %>%
    bind_rows(.id = "index") %>%
    write_tsv("./results.hla")
