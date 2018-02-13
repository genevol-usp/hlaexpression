library(Biostrings)
library(tidyverse)

sample_id <- commandArgs(TRUE)[1]
fasta <- paste0("./sample_indices/hla_", sample_id, ".fa")
index <- readDNAStringSet(fasta)

index_df <- tibble(allele = names(index), cds = as.character(index)) %>%
    filter(grepl("IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)", allele)) %>%
    mutate(start = 0L, end = nchar(cds)) %>%
    select(chrom = allele, start, end)

write_tsv(index_df, 
	  paste0("./sample_indices/", sample_id, "_imgt.bed"), 
	  col_names = FALSE)
