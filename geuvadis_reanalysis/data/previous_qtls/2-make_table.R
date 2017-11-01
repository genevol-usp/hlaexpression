devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla_v25 <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

gencode_hla_v12 <- 
    "~/gencode_data/gencode.v12.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>% 
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

gencode_hla_v19 <- 
    "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>% 
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

battle <- 
    readxl::read_excel("./battle_eqtls.xls", sheet = 1) %>%
    inner_join(gencode_hla_v25, by = c("GENE_NAME" = "gene_name")) %>%
    select(phenotype = GENE_NAME, variant = SNP_ID) %>%
    mutate(source = "Battle (2014)")

geuvadis_eur <- 
    read_tsv("./geuvadis_eur_eqtls.txt.gz", col_names = FALSE) %>%
    inner_join(gencode_hla_v12, by = c("X4" = "gene_id")) %>%
    select(phenotype = gene_name, variant = X1, slope = X10, pval = X12) %>%
    mutate(source = "Lappalainen (2013)") %>%
    left_join(rsmerge, by = c("variant" = "rsHigh")) %>%
    mutate(variant = ifelse(is.na(rsCurrent), variant, rsCurrent)) %>%
    select(source, phenotype, variant, slope, pval)

geuvadis_eur %>%
    select(-source) %>%
    write_tsv("./geuvadis_eqtl_slope_pvals.tsv")

geuvadis_eur <- select(geuvadis_eur, phenotype, variant, source)

barreiro <- 
    readxl::read_excel("./barreiro_eqtls.xlsx", 1, skip = 2) %>%
    select(phenotype = external_gene_name, variant = NI_top_snp_id) %>%
    filter(phenotype %in% gencode_hla_v25$gene_name) %>%
    mutate(source = "Nedelec (2016)")

mary <- 
    tibble(phenotype = "HLA-C", 
	   variant = c("rs2395471", "rs9264942", "rs67384697"), 
	   source = c("Vince (2016)", "Thomas (2009)", "Kulkarni (2011)"))

gtex <- 
    tibble(phenotype = c("HLA-A", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"),
	   variant = c("rs1111180", "rs1131123", "rs397709393", "rs17612907", "rs5016923"),
	   source = "GTEx")

delaneau <-
    read_delim("./LCL_eQTL.txt", col_names = FALSE, delim = " ") %>%
    inner_join(gencode_hla_v19, by = c("X1" = "gene_id")) %>%
    filter(X20 <= 0.05) %>%
    mutate(source = "Delaneau (2017)") %>%
    select(phenotype = gene_name, variant = X8, source)

rsmerge <- read_tsv("./RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, rsLow = X2, build_id = X3, rsCurrent = X7) %>%
    mutate_at(vars(rsHigh, rsLow, rsCurrent), function(x) paste0("rs", x))

qtls_df <-
    bind_rows(battle, geuvadis_eur, barreiro, mary, gtex, delaneau) %>%
    left_join(rsmerge, by = c("variant" = "rsHigh")) %>%
    mutate(variant = ifelse(is.na(rsCurrent), variant, rsCurrent)) %>%
    select(source, phenotype, variant) %>%
    arrange(source, phenotype)

write_tsv(qtls_df, "./previous_eQTL.tsv")
