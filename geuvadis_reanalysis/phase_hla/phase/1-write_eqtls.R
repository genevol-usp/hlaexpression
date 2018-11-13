library(tidyverse)

eqtls <- 
    "../eqtl_mapping/transcriptomemapping/hla_personalized/2-conditional_analysis/hla_qtls.tsv" %>%
    read_tsv() %>%
    filter(best == 1L) %>%
    pull(var_id)

write_lines(eqtls, "./eqtl_rsID.txt")
