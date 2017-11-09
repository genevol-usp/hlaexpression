devtools::load_all("~/hlaseqlib")
library(tidyverse)

extract_regdb_genes <- function(x) {

    x_split <- grep("QTL", x, value = TRUE) %>%
	strsplit(",") %>%
        unlist()

    eQTL <- x_split %>%
        grep("eQTL", ., value = TRUE) %>%
        strsplit("\\|") %>%
        map(3) %>%
        unlist()

    dsQTL <- x_split %>%
        grep("dsQTL", ., value = TRUE) %>%
        strsplit("\\|") %>%
        map(5) %>%
        unlist()

    c(eQTL, dsQTL) %>%
	unique() %>%
	.[!is.na(.)] %>%
	paste(collapse = ",")
}

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla_v25 <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>% 
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

regdb <- readLines("./RegulomeDB.dbSNP141.mhc.qtl.txt") %>%
    strsplit("\t")

regdb_df <- 
    tibble(rsid = map_chr(regdb, 4),
	   info = paste("regulomeDB:", map_chr(regdb, extract_regdb_genes)))

haploreg <- read_tsv("./eqtls_v4.1.tsv.gz", col_names = FALSE) %>%
    filter(grepl("^rs", X1)) %>%
    rename(rsid = X1, info = X2) 

battle <- readxl::read_excel("./battle_eqtls.xls", sheet = 1) %>%
    inner_join(gencode_hla_v25, by = c("GENE_NAME" = "gene_name")) %>%
    mutate(info = paste("Battle (2014):", GENE_NAME)) %>%
    select(rsid = SNP_ID, info)

barreiro <- readxl::read_excel("./barreiro_eqtls.xlsx", 1, skip = 2) %>%
    filter(external_gene_name %in% gencode_hla_v25$gene_name) %>%
    mutate(info = paste("Nedelec (2016):", external_gene_name)) %>%
    select(rsid = NI_top_snp_id, info)

mary <- 
    tibble(rsid = c("rs2395471", "rs9264942", "rs67384697"), 
	   info = paste(c("Vince (2016):", "Thomas (2009):", "Kulkarni (2011):"), 
			"HLA-C"))

delaneau <- read_delim("./LCL_eQTL.txt", col_names = FALSE, delim = " ") %>%
    inner_join(gencode_hla_v19, by = c("X1" = "gene_id")) %>%
    filter(X20 <= 0.05) %>%
    mutate(info = paste("Delaneau (2017):", gene_name)) %>%
    select(rsid = X8, info)

rsmerge <- read_tsv("./RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, rsCurrent = X7) %>%
    mutate_at(vars(rsHigh, rsCurrent), function(x) paste0("rs", x))

mhc_rsids <- readLines("./mhc_rsids_vcf.txt")

qtls_df <- bind_rows(regdb_df, haploreg, battle, barreiro, mary, delaneau) %>%
    left_join(rsmerge, by = c("rsid" = "rsHigh")) %>%
    mutate(rsid = ifelse(is.na(rsCurrent), rsid, rsCurrent)) %>%
    select(rsid, info) %>%
    filter(rsid %in% mhc_rsids)

write_tsv(qtls_df, "./qtl_catalog.tsv")
