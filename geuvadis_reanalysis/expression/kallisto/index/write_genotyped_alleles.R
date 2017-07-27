library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

index <- readDNAStringSet("/home/vitor/hlaexpression/geuvadis_reanalysis/expression/kallisto/index/imgt_index.fa")

genos <- 
  read_tsv(processed_quants) %>%
  filter(est_counts > 0) %>%
  select(subject, allele) %>%
  distinct() %>%
  separate_rows(allele, sep = "=") %>%  
  mutate(path = file.path(outdir, paste0("hla_", subject, ".fa"))) %>%
  split(.$subject)

map(genos, ~writeXStringSet(index[.$allele], unique(.$path)))
