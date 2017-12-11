devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name

samples <- geuvadis_info %>%
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications/imgt_quants.tsv")

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- imgt_quants %>%
    left_join(gencode_pri_tx, by = c("target_id" = "tx_id")) %>%
    select(tx_id = target_id, locus = gene_name, est_counts, tpm) %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize_at(vars(est_counts, tpm), sum) %>%
    ungroup()

write_tsv(out_df, "./quantifications/processed_imgt_quant.tsv") 
