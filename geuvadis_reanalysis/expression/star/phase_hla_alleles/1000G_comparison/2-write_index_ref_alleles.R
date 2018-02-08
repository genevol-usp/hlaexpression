library(tidyverse)
library(Biostrings)

index_files <- list.files(".", pattern = "index_.+\\.tsv$") 

index_df <- index_files %>%
    map_df(read_tsv) %>%
    mutate(allele = paste0("IMGT_", allele))

index <- index_df %>%
  split(.$allele) %>%
  map_chr("transcript") %>%
  DNAStringSet()

writeXStringSet(index, "./index_ref_positions.fa")
unlink(index_files)
