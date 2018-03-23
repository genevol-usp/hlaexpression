devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- readLines("../../data/sample_info/samples_phase3_ena_eur.txt")

imgt_quants <- file.path("./quantifications_winners", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name,
	   est_counts = NumReads, tpm = TPM) %>%
    arrange(subject, locus, desc(est_counts))

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd <- mutate(pag, allele = hla_trimnames(allele, 3))

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- thresholds %>%
    map_df(~hla_genotype_dt(imgt_quants, .), .id = "th")

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
    full_join(distinct(imgt_quants, locus), by = "locus") %>%
    mutate(th = ifelse(is.na(th), best_th_average, th))

geno_calls <- inner_join(typings, best_th, by = c("th", "locus")) %>%
    select(subject, locus, allele) %>%
    mutate(allele = sub("^([^=]+).*$", "\\1", allele),
	   allele = gsub("IMGT_", "", allele))

mhc_calls <- read_tsv("./quantifications_top2/genotype_calls.tsv") %>%
    mutate(allele = gsub("IMGT_", "", allele))

mhc_hom <- mhc_calls %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 1L) %>%
    ungroup()

fixed_hom <- calc_genotyping_accuracy(mhc_hom, geno_calls, by_locus = FALSE) %>%
    group_by(subject, locus) %>%
    filter(any(!correct)) %>%
    ungroup() %>%
    select(subject, locus, allele = allele.y)

fixed_loci <- distinct(fixed_hom, subject, locus)

final_genos <- anti_join(mhc_calls, fixed_loci, by = c("subject", "locus")) %>%
    bind_rows(fixed_hom) %>%
    mutate(allele = sub("^([^=]+).*$", "IMGT_\\1", allele)) %>%
    arrange(subject, locus, allele)

write_tsv(final_genos, "./quantifications_winners/genotype_calls.tsv")

final_calls <- final_genos %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("HLA-", "", locus),
	   allele = gsub("IMGT_", "", allele),
	   allele = hla_trimnames(allele, 3))

calc_genotyping_accuracy(final_calls, goldstd) %>%
    write_tsv("./genotyping_concordance.tsv")

