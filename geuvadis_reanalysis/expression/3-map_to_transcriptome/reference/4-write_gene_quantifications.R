devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

quant_dir <- "./quantifications"

gencode <- gencode_chr_tx %>%
    filter(chr %in% 1:22) %>%
    select(target_id = tx_id, gene_id, gene_name)

samples <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(name, subject = ena_id)

expression_df <- file.path(quant_dir, samples$subject, "quant.sf") %>%
    setNames(samples$subject) %>%
    map_df(read_tsv, .id = "subject") %>%
    left_join(samples, by = "subject") %>%
    select(subject = name, target_id = Name, tpm = TPM)

gene_df <- expression_df %>%
    inner_join(gencode, by = "target_id") %>%
    group_by(subject, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

write_tsv(gene_df, "gene_quantifications.tsv")
