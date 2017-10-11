devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

gencode_hla_V12 <- 
    get_gencode_coords("~/gencode_data/gencode.v12.annotation.gtf.gz") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

battle <- 
    readxl::read_excel("./battle_eqtls.xls", sheet = 1) %>%
    inner_join(gencode_hla, by = c("GENE_NAME" = "gene_name")) %>%
    select(phenotype = GENE_NAME, variant = SNP_ID) %>%
    mutate(variant = recode(variant, "rs9273448" = "rs1049225"),
	   source = "Battle (2014)")

geuvadis_eur <- 
    "../../qtls/qtls_star/imgt/rtc/geuvadis_eqtls/catalog.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    separate(X2, c("junk1", "phenotype", "junk2", "junk3"), sep = ":") %>%
    select(phenotype, variant = X1) %>%
    mutate(source = "Lappalainen (2013)")

barreiro <- 
    readxl::read_excel("./barreiro_eqtls.xlsx", 1, skip = 2) %>%
    select(phenotype = external_gene_name, variant = NI_top_snp_id) %>%
    filter(phenotype %in% gencode_hla$gene_name) %>%
    mutate(source = "Nedelec (2016)")

mary <- 
    tibble(phenotype = "HLA-C", 
	   variant = c("rs2395471", "rs9264942", "rs67384697"), 
	   source = c("Vince (2016)", "Thomas (2009)", "Kulkarni (2011)"))
 
qtls_df <-
    bind_rows(battle, geuvadis_eur, barreiro, mary) %>%
    select(source, phenotype, variant) %>%
    arrange(source, phenotype)

write_tsv(qtls_df, "./previous_eQTL.tsv")
