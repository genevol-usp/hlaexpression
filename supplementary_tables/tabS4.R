library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/functional_features/results.tsv") %>%
    write_csv("./tables/S4_table.csv")
