library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/crd/results.tsv") %>%
    select(gene, rank, eqtl = variant, crd_eqtl = crd_var, d_prime, rtc) %>%
    write_csv("./tables/S5_table.csv")
