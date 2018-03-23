library(tidyverse)

index_files <- list.files(".", pattern = "index_.+\\.tsv$") 

index_df <- map_df(index_files, read_tsv) %>%
    mutate(cds = gsub("\\*", "N", cds))

write_tsv(index_df, "./index_ref_positions.tsv")
unlink(index_files)
