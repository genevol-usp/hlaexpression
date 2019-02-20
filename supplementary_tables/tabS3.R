library(tidyverse)

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/top_previous_qtls/results.tsv") %>%
    select(gene, rank, eqtl_hlapers = qtl_personalized, 
           eqtl_previous = qtl_previous, r2, d_prime, rtc, info) %>%
    extract(info, "study", "(.+) .") %>% 
    mutate(eqtl_previous = paste0(eqtl_previous, " (", study, ")"),
           significant = as.integer(rtc > 0.9)) %>%
    select(-study) %>%
    arrange(gene, rank, desc(rtc)) %>%
    write_csv("./tables/S3_table.csv")
