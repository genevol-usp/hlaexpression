devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

bed <- select(gencode_hla, chr, start, end)

write_tsv(bed, "./hla_regions.bed", col_names = FALSE)
