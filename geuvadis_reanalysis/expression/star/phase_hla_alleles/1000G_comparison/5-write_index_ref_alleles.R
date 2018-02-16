library(tidyverse)
library(Biostrings)

index_files <- list.files(".", pattern = "index_.+\\.tsv$") 

index_df <- map_df(index_files, read_tsv)

index <- index_df %>%
    unite(id, c("allele", "exon"), sep = "_") %>% 
    filter(!is.na(cds)) %>%
    mutate(cds = gsub("\\*", "N", cds)) %>%
    split(.$id) %>%
    map_chr("cds") %>%
    DNAStringSet()

writeXStringSet(index, "./index_ref_positions.fa")
unlink(index_files)
