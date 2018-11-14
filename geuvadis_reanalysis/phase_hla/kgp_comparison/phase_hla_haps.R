devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

individual_id <- commandArgs(TRUE)[1]

hla_haps <- read_tsv("./hla_diffs_to_1000Ghaps.tsv") %>%
    filter(subject == individual_id)

phased <- phase_hla(hla_haps)

out <- paste0("./hla_haps_phased_", individual_id, ".tsv")
write_tsv(phased, out)
