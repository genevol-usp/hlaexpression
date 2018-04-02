devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- sprintf("sample_%02d", 1:50)

imgt_quants <- read_tsv("./quantifications_2/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(target_id),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = target_id, est_counts, tpm)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd_genos <- read_tsv("../../data/genos.tsv")

out_df <- hla_genotype_dt(imgt_quants, th = 0) %>%
    hla_apply_zigosity_threshold(th = 0.2)

calls <- out_df %>%
    filter(locus %in% hla_genes) %>%
    select(subject, locus, allele) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    arrange(subject, locus, allele)

accuracies <- calc_genotyping_accuracy(calls, goldstd_genos)

write_tsv(accuracies, "./genotyping_accuracies_2.tsv")

write_tsv(out_df, "./quantifications_2/processed_imgt_quants.tsv")
