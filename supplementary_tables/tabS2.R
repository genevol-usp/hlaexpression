library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/reference_qtls/results.tsv") %>%
    select(gene, rank_hlapers = rank_personalized, variant_hlapers = variant_personalized, 
	   rank_ref, variant_ref, d_prime, rtc) %>%
    mutate(d_prime = ifelse(variant_hlapers == variant_ref, "*", d_prime),
	   rtc = ifelse(variant_hlapers == variant_ref, "*", rtc)) %>%
    write_csv("./tables/S2_table.csv")
