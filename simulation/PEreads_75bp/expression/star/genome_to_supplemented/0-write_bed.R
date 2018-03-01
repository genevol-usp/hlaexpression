devtools::load_all("~/hlaseqlib")
library(tidyverse)

imgt_tx <- "~/hlaexpression/imgt_index_v2/imgt_pri_ids.txt" %>%
    readLines()

index <- 
    "~/hlaexpression/index_transcriptome/gencode.v25.PRI.uniqTranscripts.fa" %>%
    Biostrings::readDNAStringSet()

index_imgt <- index[imgt_tx]

bed <- tibble(chrom = names(index_imgt), s = as.character(index_imgt)) %>%
    mutate(start = 0L, end = nchar(s)) %>%
    select(chrom, start, end)

write_tsv(bed, "./imgt.bed", col_names = FALSE)
