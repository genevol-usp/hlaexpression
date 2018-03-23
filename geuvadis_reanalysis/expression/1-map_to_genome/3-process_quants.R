devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- readLines("../../data/sample_info/samples_phase3_ena_eur.txt")

imgt_quants <- read_tsv("./quantifications/imgt_quants.tsv")

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- imgt_quants %>%
    left_join(gencode_pri_tx, by = c("Name" = "tx_id")) %>%
    select(subject, tx_id = Name, locus = gene_name, 
	   est_counts = NumReads, tpm = TPM) %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize_at(vars(est_counts, tpm), sum) %>%
    ungroup()

write_tsv(out_df, "./quantifications/processed_imgt_quants.tsv")
