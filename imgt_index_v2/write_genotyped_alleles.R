library(Biostrings)
library(tidyverse)

processed_quants <- commandArgs(TRUE)[1]
outdir <- commandArgs(TRUE)[2]

index <- 
    readDNAStringSet("/home/vitor/hlaexpression/imgt_index_v2/imgt_index.fa")

genos <- read_tsv(processed_quants) %>%
    select(subject, allele) %>%
    mutate(path = file.path(outdir, paste0("hla_", subject, ".fa"))) %>%
    split(.$subject)

map(genos, ~writeXStringSet(index[.$allele], unique(.$path)))
