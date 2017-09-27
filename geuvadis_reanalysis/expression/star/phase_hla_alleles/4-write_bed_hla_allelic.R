devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))  

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, gene_id)

hla_by_allele <- 
    read_tsv("./data/1000G_haps_expression_snps.tsv") %>%
    filter(rank == 0L) %>%
    select(subject, locus, hap, hla_allele, tpm) 

hla_by_gene_best <-
    read_tsv("../../../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_60.bed.gz") %>%
    inner_join(gencode_hla, c("gid" = "gene_id")) %>%
    gather(subject, tpm_pc, HG00096:NA20828) %>%
    select(subject, locus = gene_name, tpm_pc) %>%
    mutate(locus = sub("^HLA-", "", locus))
  
bed_best <- 
    hla_by_allele %>%
    unite(id, locus, hap, sep = "_", remove = FALSE) %>%
    select(subject, id, locus, tpm) %>%
    group_by(subject, locus) %>%
    mutate(r = tpm/sum(tpm)) %>%
    ungroup() %>%
    left_join(hla_by_gene_best, by = c("subject", "locus")) %>%
    group_by(locus) %>%
    mutate(tpm_pc = tpm_pc + max(tpm_pc),
	   tpm = tpm_pc * r,
	   tpm = GenABEL::rntransform(tpm)) %>%
    ungroup() %>%
    select(subject, id, tpm) %>%
    spread(subject, tpm) %>%
    separate(id, c("locus", "hap"), sep = "_")

write_tsv(bed_best, "./data/hla_allele_expression_bestpc.bed")

hla_by_gene_10 <-
    read_tsv("../../../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_10.bed.gz") %>%
    inner_join(gencode_hla, c("gid" = "gene_id")) %>%
    gather(subject, tpm_pc, HG00096:NA20828) %>%
    select(subject, locus = gene_name, tpm_pc) %>%
    mutate(locus = sub("^HLA-", "", locus))
  
bed_10 <- 
    hla_by_allele %>%
    unite(id, locus, hap, sep = "_", remove = FALSE) %>%
    select(subject, id, locus, tpm) %>%
    group_by(subject, locus) %>%
    mutate(r = tpm/sum(tpm)) %>%
    ungroup() %>%
    left_join(hla_by_gene_10, by = c("subject", "locus")) %>%
    group_by(locus) %>%
    mutate(tpm_pc = tpm_pc + max(tpm_pc),
	   tpm = tpm_pc * r,
	   tpm = GenABEL::rntransform(tpm)) %>%
    ungroup() %>%
    select(subject, id, tpm) %>%
    spread(subject, tpm) %>%
    separate(id, c("locus", "hap"), sep = "_")

write_tsv(bed_10, "./data/hla_allele_expression_10pcs.bed")
