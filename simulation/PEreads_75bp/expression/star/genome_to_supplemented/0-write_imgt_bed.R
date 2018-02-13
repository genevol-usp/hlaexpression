devtools::load_all("~/hlaseqlib")
library(tidyverse)

imgt_loci <- readLines("~/hlaexpression/imgt_index_v2/imgt_loci.txt") %>%
    paste0("HLA-", .)

gencode_pri_gene %>%
    filter(gene_name %in% imgt_loci) %>%
    arrange(chr, start, end) %>%
    select(chr, start, end) %>%
    write_tsv("./hla.bed", col_names = FALSE)
