library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/functional_features/results.tsv") %>%
    mutate_at(vars(tfbs, dhs, histone_marks), ~ifelse(is.na(.), "", .)) %>%
    write_csv("./tables/S4_table.csv")
