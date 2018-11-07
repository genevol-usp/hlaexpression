library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/gwas/results.tsv") %>%
    write_csv("./tables/S6_table.csv")
