devtools::load_all("~/hlaseqlib")
library(Biostrings)
library(tidyverse)

imgt_loci <- readLines("../imgt_index/imgt_loci.txt") %>%
    paste0("HLA-", .)

mhc_tx <- gencode_pri_tx %>%
    filter(start >= 29722000L, end <= 33144000, !gene_name %in% imgt_loci) %>%
    pull(tx_id)

index <- 
    readDNAStringSet("../imgt_index/gencode.v25.PRI.IMGT.transcripts.fa") %>%
    .[grepl("IMGT", names(.)) | names(.) %in% mhc_tx]

writeXStringSet(index, "./gencode.v25.MHC.IMGT.transcripts.fa")
