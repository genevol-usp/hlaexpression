library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

index <- readDNAStringSet("/home/vitor/rnaseq/kallisto/index/imgt_index.fa")

genos <- 
  read_tsv(processed_quants) %>%
  select(subject, allele) %>%
  distinct() %>%
  separate_rows(allele, sep = "=") %>%  
  mutate(path = file.path(outdir, paste0("hla_", subject, ".fa"))) %>%
  split(.$subject)

purrr::map(genos, ~writeXStringSet(index[.$allele], unique(.$path)))
