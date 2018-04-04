devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- geuvadis_info %>% 
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

imgt_loci <- readLines("../../../../imgt_index_v2/imgt_loci.txt") %>%
    paste0("HLA-", .)

imgt_quants <- read_tsv("./quantifications/imgt_quants.tsv")

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- imgt_quants %>%
    left_join(gencode_pri_tx, by = c("Name" = "tx_id")) %>%
    select(subject, tx_id = Name, locus = gene_name, 
	   est_counts = NumReads, tpm = TPM) %>%
    filter(locus %in% imgt_loci) %>%
    group_by(subject, locus) %>%
    summarize_at(vars(est_counts, tpm), sum) %>%
    ungroup()

write_tsv(out_df, "./quantifications/processed_imgt_quants.tsv")
