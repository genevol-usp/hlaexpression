library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/crd/results.tsv") %>%
    select(gene, rank, eqtl = variant, crd_eqtl = crd_var, d_prime, rtc) %>%
    filter(rtc >= 0.95) %>%
    mutate_at(vars(d_prime, rtc), ~round(., 2)) %>%
    write_csv("./tables/S5_table.csv")
