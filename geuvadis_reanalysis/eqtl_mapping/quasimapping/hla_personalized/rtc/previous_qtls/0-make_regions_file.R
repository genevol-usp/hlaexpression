devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla %>%
    mutate(region = paste0(chr, ":", start, "-", end)) %>%
    pull(region) %>%
    write_lines("./hla_regions.txt")

