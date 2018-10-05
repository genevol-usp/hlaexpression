devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- readLines("../data/sample_ids.tsv")

imgt_quants <- read_tsv("./quantifications/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name)) %>%
    select(subject, locus, allele = Name, est_counts = NumReads, tpm = TPM)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

out_df <- hla_genotype_dt(imgt_quants, th = 0) %>%
    hla_apply_zigosity_threshold(th = 0.1)

write_tsv(out_df, "./hla_quantifications")
