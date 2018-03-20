devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- readLines("../../../data/sample_info/samples_phase3_ena_eur.txt")

quants <- file.path("./quantifications_top2", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    rename(allele = Name, est_counts = NumReads, tpm = TPM) %>%
    mutate(locus = imgt_to_gname(allele), 
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, gene_id, locus, allele, est_counts, tpm)

missing_files <- samples[! samples %in% quants$subject]

if (length(missing_files) > 0L) {
     stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd <- mutate(pag, allele = hla_trimnames(allele, 3))

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- thresholds %>%
    map_df(~hla_genotype_dt(quants, .), .id = "th")

calls <- typings %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(th, subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(as.character(subject)),
           locus = sub("^HLA-", "", locus),
           allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    arrange(subject, locus, allele)

accuracies <- calls %>%
    split(.$th) %>%
    map_df(~calc_genotyping_accuracy(., goldstd), .id = "th") %>%
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
    full_join(distinct(quants, locus), by = "locus") %>%
    mutate(th = ifelse(is.na(th), best_th_average, th))

geno_calls <- inner_join(typings, best_th, by = c("th", "locus")) %>%
    select(subject, locus, allele)

write_tsv(geno_calls, "./quantifications_top2/genotype_calls.tsv")
