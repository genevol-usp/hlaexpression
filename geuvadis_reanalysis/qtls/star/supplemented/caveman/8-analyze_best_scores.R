devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- select(gencode_hla, gene_id, gene_name)

imgt_qtls <- 
    read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    #group_by(phen_id, var_id) %>%
    #filter(bwd_pval == min(bwd_pval)) %>% 
    #group_by(phen_id, var_id) %>%
    #filter(rank == min(rank)) %>%
    #ungroup() %>%
    select(gene_id = phen_id, gene_name, rank, var_id, pos = var_from)

results_all_imgt <- data.table::fread("./results.all") %>% as_tibble()

max_scores <- results_all_imgt %>% 
    group_by(GENE) %>%
    filter(CaVEMaN == max(CaVEMaN)) %>%
    ungroup() %>%  
    separate(GENE, c("gene_id", "eqtl_chrom", "eqtl_pos", "ref", "alt"), sep = "_") %>%
    filter(gene_id %in% imgt_qtls$gene_id, eqtl_pos != POS) %>%
    left_join(imgt_qtls, by = c("gene_id", "POS" = "pos")) %>%
    select(gene_name, rank, var_id, eqtl_pos, caveman_var_pos = POS, cor = COR, 
	   pvalue = P, score = CaVEMaN)

