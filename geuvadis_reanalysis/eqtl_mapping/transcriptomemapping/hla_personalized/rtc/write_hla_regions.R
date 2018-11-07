devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

gencode_hla %>%
    select(chr, start, end) %>%
    unite(coord, start:end, sep = "-") %>%
    unite(info, chr:coord, sep = ":") %>%
    write_lines("./hla_regions.txt")
