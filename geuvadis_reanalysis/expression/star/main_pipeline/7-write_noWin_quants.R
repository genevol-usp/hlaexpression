library(tidyverse)

samples <- readLines("../../../data/sample_info/samples_phase3_ena_eur.txt")

quants <- file.path("./quantifications_top5", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    rename(allele = Name, est_counts = NumReads) %>%
    mutate(locus = sub("^IMGT_([^\\*]+).+$", "HLA-\\1", allele)) %>%
    select(subject, locus, allele, est_counts)

winners <- quants %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

winners_list <- winners %>%
    select(subject, allele) %>%
    mutate(path = file.path("./quantifications_top5", subject, "winner_alleles.txt")) %>%
    split(.$subject)

map(winners_list, ~writeLines(.$allele, unique(.$path)))
