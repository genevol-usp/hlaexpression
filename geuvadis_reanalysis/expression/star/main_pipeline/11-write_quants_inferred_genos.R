library(tidyverse)

samples <- readLines("../../../data/sample_info/samples_phase3_ena_eur.txt")

quants <- file.path("./quantifications_noWin", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    rename(allele = Name, est_counts = NumReads) %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "HLA-\\1", allele)) %>%
    select(subject, locus, allele, est_counts)

winners_2nd <- quants %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

winners <- read_tsv("./quantifications_top5/winners_quants.tsv")

bind_rows(winners, winners_2nd) %>%
    arrange(subject, locus, desc(est_counts)) %>%
    write_tsv("./quants_inferred_genotypes.tsv")
