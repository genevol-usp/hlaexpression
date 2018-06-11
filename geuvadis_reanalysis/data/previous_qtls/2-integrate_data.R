devtools::load_all("~/hlaseqlib")
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
    select(gene_name, rsid = SNP_ID, pval = log10pvalue, min_rank) %>%
    arrange(gene_name, min_rank)

geuvadis_exon <- read_tsv("./EUR373.exon.cis.FDR5.all.rs137.txt.gz") %>%
    inner_join(gencode_hla_v12, by = c("GENE_ID" = "gene_id")) %>%
    group_by(GENE_ID) %>%
    mutate(min_rank = rank(pvalue, ties = "min")) %>%
    ungroup() %>%
    select(gene_name, rsid = SNP_ID, pval = log10pvalue, min_rank) %>%
    arrange(gene_name, min_rank)

gtex_gene <-
    paste0("./gtex_", sub("HLA-", "", gencode_hla$gene_name), "_V7.csv") %>%
    map_df(. %>% read_csv() %>%
	   filter(Tissue == "Cells - EBV-transformed lymphocytes", 
		  `Gencode Id` %in% gencode_hla_v19$gene_id) %>%
	   mutate(min_rank = rank(`P-Value`, ties = "min"),
		  pval = -log10(`P-Value`)) %>%
	   select(gene_name = `Gene Symbol`, rsid = `SNP Id`, pval, min_rank)) %>%
    arrange(gene_name, min_rank)
    
delaneau2018 <- read_delim("./delaneau_LCL_eqtls.txt", col_names = FALSE, delim = " ") %>%
    inner_join(gencode_hla_v19, by = c("X1" = "gene_id")) %>%
    filter(X20 <= 0.05) %>%
    group_by(X1) %>%
    mutate(min_rank = rank(X20, ties = "min")) %>%
    ungroup() %>%
    mutate(pval = -log10(X20)) %>%
    select(gene_name, rsid = X8, pval, min_rank)

vince2017 <- tibble(gene_name = "HLA-C",
		rsid = "rs2395471",
		pval = NA,
		min_rank = NA)

thomas2009 <- tibble(gene_name = "HLA-C",
		     rsid = "rs9264942",
		     pval = NA,
		     min_rank = NA)

kulkarni2011 <- tibble(gene_name = "HLA-C",
		       rsid = "rs67384697",
		       pval = NA,
		       min_rank = NA)

petersdorf2015 <- tibble(gene_name = "HLA-DPB1",
			 rsid = "rs9277534",
			 pval = NA,
			 min_rank = NA)

integrated_data <-
    list(geuvadis_gene = geuvadis_gene,
	 geuvadis_exon = geuvadis_exon,
	 gtex_v7 = gtex_gene,
	 delaneau2018 = delaneau2018,
	 thomas2009 = thomas2009,
	 kulkarni2011 = kulkarni2011,
	 vince2017 = vince2017,
	 petersdorf2015 = petersdorf2015) %>%
    bind_rows(.id = "study") %>%
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
