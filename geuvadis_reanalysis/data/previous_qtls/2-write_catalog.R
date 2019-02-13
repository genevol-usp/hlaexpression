devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

convert_rsIDs <- function(df) {
    left_join(df, rsmerge, by = c("rsid" = "rsHigh")) %>%
    mutate(rsid = ifelse(is.na(rsCurrent), rsid, rsCurrent)) %>%
    select(names(df))
}

rsmerge <- read_tsv("./RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, build_id = X3, rsCurrent = X7) %>%
    filter(build_id <= 149) %>%
    mutate_at(vars(rsHigh, rsCurrent), function(x) paste0("rs", x)) %>%
    select(rsHigh, rsCurrent)

gencode_hla_v12 <- select(gencode_hla_v12, gene_id, gene_name)
gencode_hla_v19 <- select(gencode_hla_v19, gene_id, gene_name)

geuvadis_gene <- read_tsv("./EUR373.gene.cis.FDR5.all.rs137.txt.gz") %>%
    inner_join(gencode_hla_v12, by = c("GENE_ID" = "gene_id")) %>%
    group_by(GENE_ID) %>%
    mutate(min_rank = rank(pvalue, ties = "min")) %>%
    ungroup() %>%
    mutate(study = "GeuvadisGene") %>%
    select(study, gene_name, rsid = SNP_ID, pval = log10pvalue, min_rank) %>%
    arrange(gene_name, min_rank)

geuvadis_exon <- read_tsv("./EUR373.exon.cis.FDR5.all.rs137.txt.gz") %>%
    inner_join(gencode_hla_v12, by = c("GENE_ID" = "gene_id")) %>%
    group_by(GENE_ID) %>%
    mutate(min_rank = rank(pvalue, ties = "min")) %>%
    ungroup() %>%
    mutate(study = "GeuvadisExon") %>%
    select(study, gene_name, rsid = SNP_ID, pval = log10pvalue, min_rank) %>%
    arrange(gene_name, min_rank)

gtex_gene <-
    paste0("./gtex_", sub("HLA-", "", gencode_hla$gene_name), "_V7.csv") %>%
    map_df(. %>% read_csv() %>%
	   filter(Tissue == "Cells - EBV-transformed lymphocytes", 
		  `Gencode Id` %in% gencode_hla_v19$gene_id) %>%
	   mutate(min_rank = rank(`P-Value`, ties = "min"),
		  pval = -log10(`P-Value`)) %>%
	   select(gene_name = `Gene Symbol`, rsid = `SNP Id`, pval, min_rank)) %>%
    mutate(study = "GTExV7") %>%
    select(study, everything()) %>%
    arrange(gene_name, min_rank)
    
delaneau2017 <- read_delim("./delaneau_LCL_eqtls.txt", col_names = FALSE, delim = " ") %>%
    inner_join(gencode_hla_v19, by = c("X1" = "gene_id")) %>%
    filter(X20 <= 0.05) %>%
    group_by(X1) %>%
    mutate(min_rank = rank(X20, ties = "min")) %>%
    ungroup() %>%
    mutate(study = "Delaneau2017",
	   pval = -log10(X20)) %>%
    select(study, gene_name, rsid = X8, pval, min_rank)

fairfax2012 <- readxl::read_excel("./ng.2205-S2.xls", 3) %>%
    filter(Gene %in% gencode_hla_v19$gene_name) %>%
    select(gene_name = Gene, rsid = SNP, pval = P) %>%
    group_by(gene_name) %>%
    mutate(min_rank = rank(pval, ties = "min")) %>%
    ungroup() %>%
    bind_rows(tibble(gene_name = "HLA-C", rsid = "rs10484554", 
		     pval = 5.2e-16, min_rank = 1)) %>%
    mutate(study = "Fairfax2012") %>%
    select(study, everything()) %>%
    arrange(gene_name, min_rank)


vince2016 <- tibble(study = "Vince2016",
		    gene_name = "HLA-C",
		    rsid = "rs2395471",
		    pval = NA,
		    min_rank = NA)

thomas2009 <- tibble(study = "Thomas2009",
		     gene_name = "HLA-C",
		     rsid = "rs9264942",
		     pval = NA,
		     min_rank = NA)

thomas2012 <- tibble(study = "Thomas2012",
		     gene_name = "HLA-DPB1",
		     rsid = "rs9277534",
		     pval = NA,
		     min_rank = NA)

kulkarni2011 <- tibble(study = "Kulkarni2011",
		       gene_name = "HLA-C",
		       rsid = "rs67384697",
		       pval = NA,
		       min_rank = NA)

raj2016 <- tibble(study = "Raj2016",
		  gene_name = c("HLA-DQA1", "HLA-DQB1", "HLA-DRB1"),
		  rsid = "rs9271593",
		  pval = NA,
		  min_rank = NA)

ou2019 <- tibble(study = "Ou2019",
		 gene_name = "HLA-DPA1",
		 rsid = "rs3077",
		 pval = NA,
		 min_rank = NA)

integrated_data <-
    bind_rows(geuvadis_gene, geuvadis_exon, gtex_gene, delaneau2017, 
	      fairfax2012, thomas2009, thomas2012, kulkarni2011, vince2016, 
	      raj2016, ou2019) %>%
    convert_rsIDs()

eqtl_catalog <- integrated_data %>%
    group_by(gene_name, rsid, study) %>%
    filter(pval == max(pval)) %>%
    ungroup() %>%
    unite(info, c("study", "gene_name", "pval"), sep = "|") %>%
    group_by(rsid) %>%
    summarise(info = paste(info, collapse = ";")) %>%
    ungroup()

write_tsv(eqtl_catalog, "./qtl_catalog.tsv", col_names = FALSE)

top_eqtl_catalog <- integrated_data %>%
    group_by(gene_name, rsid, study) %>%
    filter(min_rank == 1L | is.na(min_rank)) %>%
    ungroup() %>%
    unite(info, c("study", "gene_name", "pval"), sep = "|") %>%
    group_by(rsid) %>%
    summarise(info = paste(info, collapse = ";")) %>%
    ungroup()

write_tsv(top_eqtl_catalog, "./top_qtl_catalog.tsv", col_names = FALSE)
