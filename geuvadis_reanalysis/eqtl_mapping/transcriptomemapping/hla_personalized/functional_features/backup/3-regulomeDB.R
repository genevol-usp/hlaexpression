library(tidyverse)
library(haploR)
		      
qtls <- read_tsv("../2-conditional_analysis/hla_qtls.tsv") %>%
    filter(best == 1) %>%
    select(gene, rsid = var_id, rank)

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

