library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/reference_qtls/results.tsv") %>%
    select(gene, rank_hlapers = rank_personalized, eqtl_hlapers = variant_personalized, 
	   rank_ref, eqtl_ref = variant_ref, d_prime, rtc) %>%
    mutate(d_prime = round(d_prime, 2),
           d_prime = ifelse(eqtl_hlapers == eqtl_ref, "*", d_prime),
           rtc = ifelse(eqtl_hlapers == eqtl_ref, "*", rtc)) %>%
    write_csv("./tables/S2_table.csv")
