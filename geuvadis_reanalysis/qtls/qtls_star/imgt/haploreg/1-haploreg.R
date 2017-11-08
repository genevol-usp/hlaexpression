devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(haploR)
		      
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

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

qtls <-
    read_qtltools("../3-conditional_analysis/conditional_60_all.txt.gz") %>%
    filter(bwd_best == 1) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rsID = var_id, rank)

haploreg_query <-
    queryHaploreg(query = qtls$variant, ldThresh = 1, ldPop = "EUR", 
		  epi = "vanilla", cons = "siphy", genetypes = "gencode") %>%
    filter(rsID == query_snp_rsid) %>%
    select(rsID, pos = pos_hg38, Chromatin_States, Chromatin_Marks, DNAse, 
	   Proteins, eQTL, gwas, grasp, Motifs, GENCODE_name, GENCODE_distance, 
	   dbSNP_functional_annotation, Promoter_histone_marks, 
	   Enhancer_histone_marks)

same_gene_df <- 
    haploreg_query %>%
    left_join(qtls, by = "rsID") %>%
    select(rsID, gene, rank, eQTL) %>% 
    mutate(eQTL = map(eQTL, eQTL_to_df)) %>%
    unnest() %>%
    filter(gene == gene1) %>%
    select(-gene1) %>%
    left_join(qtls, ., by = c("gene", "rank", "rsID")) %>%
    arrange(gene, rank, desc(pvalue))

write_tsv(same_gene_df, "./haploreg_results.tsv")
