library(tidyverse)
devtools::load_all("~/hlaseqlib")

gencode_hla_tx <- filter(gencode_pri_tx, gene_name %in% gencode_hla$gene_name) %>%
    select(tx_id, gene_name)

phenotypes <- read_tsv("./data/phenotypes_counts_tpm.tsv") %>%
    inner_join(gencode_hla_tx, by = c("Name" = "tx_id")) %>%
    select(gene = gene_name, tx_id = Name, TrueTPM)

phenotypes_by_gene <- phenotypes %>%
    group_by(gene) %>%
    summarise(`tpm (true)` = sum(TrueTPM)) %>%
    ungroup()

est_cds_index <- 
    read_tsv("./expression/star/quantifications_2/imgt_quants.tsv") %>%
    mutate(gene_name = imgt_to_gname(Name)) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene = gene_name, `tpm (cds)` = TPM)

est_cds_utr_index <-
    read_tsv("./expression/star_utr_index/quantifications_2/imgt_quants.tsv") %>%
    mutate(gene_name = imgt_to_gname(Name)) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene = gene_name, `tpm (cds+utr)` = TPM)

est_trimmed_index <-
    read_tsv("./expression/star_trimmed_index/quantifications_2/imgt_quants.tsv") %>%
    mutate(gene_name = imgt_to_gname(Name)) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene = gene_name, `tpm (trimmed)` = TPM)

out <- 
    left_join(phenotypes_by_gene, est_cds_index, by = "gene") %>%
    left_join(est_cds_utr_index, by = "gene") %>%
    left_join(est_trimmed_index, by = "gene") %>%
    mutate(`estimated/true (cds)` = `tpm (cds)`/`tpm (true)`,
	   `estimated/true (cds+utr)` = `tpm (cds+utr)`/`tpm (true)`,
	   `estimated/true (trimmed)` = `tpm (trimmed)`/`tpm (true)`)

write_tsv(out, "./results.tsv")
