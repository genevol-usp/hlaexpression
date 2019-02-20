library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/gwas/results.tsv") %>%
    mutate(`trait (GWAS variant)` = ifelse(`trait (GWAS variant)` == "NA (NA)", "",  `trait (GWAS variant)`)) %>%
    write_csv("./tables/S6_table.csv")
