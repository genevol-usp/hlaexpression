library(tidyverse)

ids <- 
    "/home/vitor/hlaexpression/geuvadis_reanalysis/data/sample_info/samples_phase3_eur.txt" %>%
    readLines()

result_files <- paste0("hla_haps_phased_", ids, ".tsv") %>%
    setNames(ids)

phase_df <- map_df(result_files, read_tsv, .id = "subject")

unlink(result_files)
write_tsv(phase_df, "hla_haps_phased.tsv")
