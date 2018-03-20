devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- geuvadis_info %>%
    filter(kgp_phase3 == 1L, pop != "YRI") %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications_MHC/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name, 
	   est_counts = NumReads, tpm = TPM)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd_genos <- mutate(pag, allele = hla_trimnames(allele, 3))

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- thresholds %>%
    map_df(~hla_genotype_dt(imgt_quants, .), .id = "th")

calls <- typings %>%
    filter(locus %in% hla_genes) %>%
    select(th, subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("HLA-", "", locus),
	   allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    arrange(subject, locus, allele)

accuracies <- calls %>%
    split(.$th) %>%
    map_df(~calc_genotyping_accuracy(., goldstd_genos), .id = "th") %>%
    group_by(th) %>%
    mutate(th_average = mean(accuracy)) %>%
    ungroup()

best_th_average <- accuracies %>%
  slice(which.max(th_average)) %>%
  pull(th) %>%
  as.character()

best_th <- accuracies %>%
    group_by(locus) %>%
    slice(which.max(accuracy)) %>%
    ungroup() %>%
    mutate(locus = paste0("HLA-", locus),
	   th = as.character(th)) %>%
    select(th, locus) %>%
    full_join(distinct(imgt_quants, locus), by = "locus") %>%
    mutate(th = ifelse(is.na(th), best_th_average, th))

out_df <- inner_join(typings, best_th, by = c("th", "locus")) %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus, allele)

write_tsv(out_df, "./quantifications_MHC/allele_calls.tsv")
