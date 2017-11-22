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
    
convert_rsIDs <- function(df) {
    
    left_join(df, rsmerge, by = c("rsid" = "rsHigh")) %>%
    mutate(rsid = ifelse(is.na(rsCurrent), rsid, rsCurrent)) %>%
    select(rsid, info)
}

hla_genes <- sort(gencode_hla$gene_name)

mhc_rsids <- readLines("./mhc_rsids_vcf.txt")

gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>% 
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

rsmerge <- read_tsv("./RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, rsCurrent = X7) %>%
    mutate_at(vars(rsHigh, rsCurrent), function(x) paste0("rs", x))

# regulomeDB
regdb <- readLines("./RegulomeDB.dbSNP141.mhc.qtl.txt") %>%
    strsplit("\t")

regdb_df <- 
    tibble(rsid = map_chr(regdb, 4),
	   info = paste("regulomeDB:", map_chr(regdb, extract_regdb_genes))) %>%
    convert_rsIDs() %>%
    filter(rsid %in% mhc_rsids) 

# haploreg
# read with fread bc the file as one \r within a line and somehow fread handles
# this
haploreg <- 
    data.table::fread("zcat < ./eqtls_v4.1.tsv.gz", header = FALSE, sep = "\t") %>%
    as_tibble() %>%
    filter(V1 != ".") %>%
    mutate(V1 = ifelse(grepl("AFFX", V1), sub("^.+(rs\\d+)$", "\\1", V1), V1)) %>%
    rename(rsid = V1, info = V2) %>% 
    convert_rsIDs() %>%
    filter(rsid %in% mhc_rsids) %>%
    separate_rows(info, sep = ";") 
    
haploreg_comma <- filter(haploreg, grepl("\"", info)) %>%
    extract(info, c("study", "tissue", "gene", "pvalue"), 
	    "([^,]+),([^,]+),\"([^\"]+)\",([^,]+)") %>%
    mutate(gene = gsub(",", "/", gene))

haploreg_v <- filter(haploreg, !grepl("\"", info)) %>% 
    separate(info, c("study", "tissue", "gene", "pvalue"), sep = ",")

haploreg <- bind_rows(haploreg_v, haploreg_comma) %>%
    mutate(pvalue = sub("\\\r$", "", pvalue),
	   pvalue = as.numeric(pvalue)) %>%
    group_by(study, tissue, gene) %>%
    filter(pvalue == min(pvalue)) %>%
    ungroup() %>%
    mutate(pvalue = trimws(format(pvalue, scientific = TRUE))) %>%
    unite(info, study:pvalue, sep = ",") %>%
    group_by(rsid) %>%
    summarize(info = paste(info, collapse = ";")) %>%
    ungroup()

# selected pubs:
battle <- readxl::read_excel("./battle_eqtls.xls", sheet = 1) %>%
    inner_join(gencode_hla, by = c("GENE_NAME" = "gene_name")) %>%
    mutate(info = paste("Battle (2014):", GENE_NAME)) %>%
    select(rsid = SNP_ID, info) %>%
    convert_rsIDs() %>%
    filter(rsid %in% mhc_rsids)

barreiro <- readxl::read_excel("./barreiro_eqtls.xlsx", 1, skip = 2) %>%
    filter(external_gene_name %in% gencode_hla$gene_name) %>%
    mutate(info = paste("Nedelec (2016):", external_gene_name)) %>%
    select(rsid = NI_top_snp_id, info) %>%
    convert_rsIDs() %>%
    filter(rsid %in% mhc_rsids)

#mary <- 
#    tibble(rsid = c("rs2395471", "rs9264942", "rs67384697"), 
#	   info = paste(c("Vince (2016):", "Thomas (2009):", "Kulkarni (2011):"), 
#			"HLA-C"))

delaneau <- read_delim("./LCL_eQTL.txt", col_names = FALSE, delim = " ") %>%
    inner_join(gencode_hla_v19, by = c("X1" = "gene_id")) %>%
    filter(X20 <= 0.05) %>%
    mutate(info = paste("Delaneau (2017):", gene_name)) %>%
    select(rsid = X8, info) %>%
    convert_rsIDs() %>%
    filter(rsid %in% mhc_rsids)

# merge:
qtls_df <- bind_rows(regdb_df, haploreg, battle, barreiro, delaneau) %>%
    distinct(rsid, info) %>%
    group_by(rsid) %>%
    summarize(info = paste(info, collapse = ";")) %>%
    ungroup()

write_tsv(qtls_df, "./qtl_catalog.tsv")

