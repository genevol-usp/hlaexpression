devtools::load_all("~/genomicRutils")
library(tidyverse)

geuvadis_samples_phase3 <-
  geuvadis_info %>%
  filter(kgp_phase3 == 1L) %>%
  select(subject = name, ena_id, pop) %>%
  arrange(ena_id)

geuvadis_samples_phase3$ena_id %>%
  writeLines("samples_phase3_ena.txt")

geuvadis_samples_phase3$subject %>%
  sort() %>%
  writeLines("samples_phase3.txt")
