library(tidyverse)

hla_diffs <- read_tsv("./hla_diffs_to_1000Ghaps.tsv")

hla_haps <- read_tsv("./hla_haps_phased.tsv")

confi_set <- hla_haps %>%
    left_join(hla_diffs) %>%
    group_by(subject, locus) %>%
    filter(all(diffs <= 6)) %>%
    ungroup()

write_tsv(confi_set, "./hla_haps_phased_confidence_set.tsv") 
