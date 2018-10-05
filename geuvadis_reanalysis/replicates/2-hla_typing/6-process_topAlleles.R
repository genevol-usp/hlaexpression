devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

samples <- readLines("../data/sample_ids.tsv")

quants <- file.path("./quantifications_topAlleles", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    mutate(locus = imgt_to_gname(Name),
	   allele = sub("IMGT_", "", Name),
	   lineage = hla_trimnames(allele, 1)) %>%
    select(subject, locus, allele, lineage, est_counts = NumReads, tpm = TPM)

lineage_typing <- quants %>%
    group_by(subject, locus, allele = lineage) %>%
    summarise_at(vars(est_counts, tpm), sum) %>%
    ungroup() %>%
    hla_genotype_dt(th = .15) %>%
    select(subject, locus, lineage = allele) %>%
    distinct()

quants2 <- inner_join(quants, lineage_typing)

hets <- quants2 %>%
    group_by(subject, locus) %>%
    filter(n_distinct(lineage) > 1L) %>%
    group_by(subject, lineage) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

quants3 <- quants2 %>%
    group_by(subject, locus) %>%
    filter(n_distinct(lineage) == 1L) %>%
    ungroup() %>%
    bind_rows(hets) %>%
    arrange(subject, locus)

allele_typing <- hla_genotype_dt(quants3, th = .15)

sample_info <- read_tsv("../data/sample_info.tsv")

calls <- allele_typing %>%
    left_join(sample_info, by = c("subject" = "ena_id")) %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    select(subject = subject_id, locus, allele)

goldstd <- mutate(pag, allele = hla_trimnames(allele, 3))

concordance <- calc_genotyping_accuracy(calls, goldstd)

allele_typing %>%
    select(subject, locus, allele) %>%
    mutate(allele = sub("^([^=]+).*$", "\\1", allele)) %>%
    left_join(quants3) %>%
    select(subject, locus, allele, est_counts, tpm) %>%
    write_tsv("./quantifications_topAlleles/processed_imgt_quants.tsv")

write_tsv(concordance, "./quantifications_topAlleles/genotyping_concordance.tsv")
