devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

geuvadis_samples_phase3 <-
  geuvadis_info %>%
  filter(kgp_phase3 == 1L) %>%
  select(name, ena_id, pop)

geuvadis_samples_phase3 %>%
  pull(name) %>%
  writeLines("samples_phase3.txt")

geuvadis_samples_phase3 %>%
  pull(ena_id) %>%
  sort() %>%
  writeLines("samples_phase3_ena.txt")

geuvadis_samples_phase3 %>%
  filter(pop != "YRI") %>%
  pull(ena_id) %>%
  sort() %>%
  writeLines("samples_phase3_ena_eur.txt")

