devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- readLines("../data/sample_ids.tsv")

imgt_quants <- read_tsv("./quantifications_MHC/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name)) %>%
    select(subject, locus, allele = Name, est_counts = NumReads, tpm = TPM)

top_alleles <- imgt_quants %>%
    group_by(subject, locus) %>%
    top_n(5, est_counts) %>%
    ungroup() %>%
    mutate(lineage = hla_trimnames(sub("IMGT_", "", allele), 1)) %>%
    group_by(subject, locus, lineage) %>%
    filter(tpm/max(tpm) > 0.25) %>%
    ungroup() %>%
    select(subject, locus, allele, est_counts, tpm)

write_tsv(top_alleles, "./quantifications_MHC/imgt_quants_topAlleles.tsv")
