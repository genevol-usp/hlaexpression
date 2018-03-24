devtools::load_all("/home/vitor/genomicRutils")

geuvadis_ph3 <- 
  dplyr::filter(geuvadis_info, kgp_phase3 == 1) %>%
  dplyr::select(name, pop)

geuvadis_ph3 %>%
  dplyr::pull(name) %>%
  writeLines("samples.all")

geuvadis_ph3 %>%
  dplyr::filter(pop != "YRI" ) %>%
  dplyr::pull(name) %>%
  writeLines("samples.eur")

geuvadis_ph3 %>%
  dplyr::filter(pop == "YRI" ) %>%
  dplyr::pull(name) %>%
  writeLines("samples.yri")
