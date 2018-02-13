devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_haps <- read_tsv("./hla_haps_mapped_to_1000G.tsv")

phased <- plyr::ddply(hla_haps, ~subject, . %>% phase_hla, .id = "subject")

write_tsv(phased, "./hla_haps_mapped_to_1000G_phased.tsv")
