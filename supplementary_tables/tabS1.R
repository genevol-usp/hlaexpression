library(tidyverse)

star <- read_tsv("../geuvadis_reanalysis/expression/2-hla_typing/genotyping_concordance.tsv") %>%
    rename(`STAR-Salmon` = accuracy)

kallisto <- read_tsv("../geuvadis_reanalysis/expression/5-pseudoalignment/hla_personalized/genotyping_accuracies_2.tsv") %>%
    rename(kallisto = accuracy)

left_join(star, kallisto) %>%
    mutate(locus = paste0("HLA-", locus)) %>%
    mutate_at(vars(2:3), ~paste0(round(.*100, 1), "%")) %>%
    write_csv("./tables/S1_table.csv")
