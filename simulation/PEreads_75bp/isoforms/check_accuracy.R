library(tidyverse)
devtools::load_all("~/hlaseqlib")

gencode_hla_tx <- filter(gencode_pri_tx, gene_name %in% gencode_hla$gene_name) %>%
    select(tx_id, gene_name)

phenotypes <- read_tsv("./data/phenotypes_counts_tpm.tsv") %>%
    inner_join(gencode_hla_tx, by = c("Name" = "tx_id")) %>%
    select(gene_name, tx_id = Name, TrueTPM)

phenotypes_by_gene <- phenotypes %>%
    group_by(gene_name) %>%
    summarise(true_tpm = sum(TrueTPM)) %>%
    ungroup()

est_star <- 
    read_tsv("./expression/star/quantifications_2/imgt_quants.tsv") %>%
    mutate(gene_name = imgt_to_gname(Name)) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, est_tpm_star = TPM)

est_star_utr_index <-
    read_tsv("./expression/star_utr_index/quantifications_2/imgt_quants.tsv") %>%
    mutate(gene_name = imgt_to_gname(Name)) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, est_tpm_star_utr_index = TPM)

out <- left_join(phenotypes_by_gene, est_star, by = "gene_name") %>%
    left_join(est_star_utr_index, by = "gene_name") %>%
    mutate(proportion_true_est_star = est_tpm_star/true_tpm,
	   proportion_true_est_star_utr_index = est_tpm_star_utr_index/true_tpm)

write_tsv(out, "./results.tsv")
