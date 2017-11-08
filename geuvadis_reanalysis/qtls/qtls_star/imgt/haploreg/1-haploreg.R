devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(haploR)
		      
hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

qtls <-
    read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rsid = var_id, rank)

# haploreg
eQTL_to_df <- function(x) {
    as.character(x) %>%
    strsplit(";") %>%
    map(as.data.frame) %>%
    bind_rows() %>%
    setNames("X1") %>% 
    separate(X1, c("study", "tissue", "gene", "pvalue"), sep = ",") %>%
    mutate(gene = ifelse(grepl("^ENSG", gene), 
			 sub("^([^_.]+).*$", "\\1", gene), 
			 gene),
	   gene = ifelse(grepl("^ENSG", gene),
			 gencode_noversion$gene_name[match(gene, gencode_noversion$gene_id)],
			 gene),
	   pvalue = -log10(as.numeric(pvalue))) %>%
    group_by(study, tissue, gene) %>%
    filter(pvalue == max(pvalue)) %>%
    ungroup()
}

gencode_noversion <- gencode_chr_gene %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id))

haploreg_query <-
    queryHaploreg(query = qtls$rsid, ldThresh = 1, ldPop = "EUR", 
		  epi = "vanilla", cons = "siphy", genetypes = "gencode") %>%
    filter(rsID == query_snp_rsid) %>%
    select(rsid = rsID, pos = pos_hg38, Chromatin_States, Chromatin_Marks, 
	   DNAse, Proteins, eQTL, gwas, grasp, Motifs, GENCODE_name, 
	   GENCODE_distance, dbSNP_functional_annotation, 
	   Promoter_histone_marks, Enhancer_histone_marks)

same_gene_df <- 
    haploreg_query %>%
    left_join(qtls, by = "rsid") %>%
    select(rsid, gene, rank, eQTL) %>% 
    mutate(eQTL = map(eQTL, eQTL_to_df)) %>%
    unnest() %>%
    filter(gene == gene1) %>%
    select(-gene1) %>%
    left_join(qtls, ., by = c("gene", "rank", "rsid")) %>%
    arrange(gene, rank, desc(pvalue))

write_tsv(same_gene_df, "./haploreg_results.tsv")

# regulomeDB
regDB <- queryRegulome(query = qtls$rsid, check_bad_snps = TRUE)

regDB_hits <- 
    regDB$res.table %>%
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
    drop_na(qtl) %>%
    select(rsid, score, qtl) %>%
    left_join(qtls, ., by = "rsid") %>%
    arrange(gene, rank) %>%
    write_tsv("./regulomeDB_results.tsv")


