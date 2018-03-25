devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(haploR)
		      
hla_genes <- gencode_hla$gene_name

gencode_hla <- select(gencode_hla, gene_id, gene_name)

qtls <-
    read_qtltools("../3-conditional_analysis/conditional_50_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rsid = var_id, rank)

regDB <- queryRegulome(query = qtls$rsid, check_bad_snps = TRUE)

regDB_hits <- regDB$res.table %>%
    select(rsid, score, hits) %>%
    separate_rows(hits, sep = ",") %>%
    mutate(hits = trimws(hits)) %>%
    extract(hits, c("phenotype", "info"), "([^\\|]+)\\|(.*)") %>%
    mutate(info = gsub("^\\||\\|$", "", info)) %>%
    group_by(rsid, score, phenotype) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    spread(phenotype, info) %>%
    select(rsid, score, chromatin_struct = Chromatin_Structure,
	   motifs = Motifs, protein_binding = Protein_Binding, 
	   qtl = Single_Nucleotides)

regDB_hits %>%
    group_by(rsid, score) %>%
    summarize(qtl = paste(unique(qtl), collapse = ";")) %>%
    ungroup() %>%
    left_join(qtls, ., by = "rsid") %>%
    arrange(gene, rank) %>%
    write_tsv("./regulomeDB_results.tsv")
