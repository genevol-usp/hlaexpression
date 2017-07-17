devtools::load_all("/home/vitor/genomicRutils")
library(tidyverse)

samples <- geuvadis_info %>%
  filter(kgp_phase3 == 1L) %>%
  pull(name)

labs <- levels(reorder(geuvadis_info$lab, geuvadis_info$lab_code))
pops <- unique(geuvadis_info$pop)

matrix_lab <- 
  matrix(0, nrow = length(samples), ncol = length(labs), 
	 dimnames = list(NULL, labs)) %>%
  as_tibble()

matrix_pop <- 
  matrix(0, nrow = length(samples), ncol = length(pops), 
	 dimnames = list(NULL, pops)) %>%
  as_tibble()

covs <- 
  geuvadis_info %>%
  filter(name %in% samples) %>%
  select(name, lab, pop) %>%
  bind_cols(matrix_lab, matrix_pop) %>%
  mutate(UNIGE = as.integer(lab == "UNIGE"),
	 CNAG_CRG = as.integer(lab == "CNAG_CRG"),
	 MPIMG = as.integer(lab == "MPIMG"),
	 ICMB = as.integer(lab == "ICMB"),
	 HMGU = as.integer(lab == "HMGU"),
	 UU = as.integer(lab == "UU"),
	 LUMC = as.integer(lab == "LUMC"),
	 GBR = as.integer(pop == "GBR"),
	 FIN = as.integer(pop == "FIN"),
	 CEU = as.integer(pop == "CEU"),
	 YRI = as.integer(pop == "YRI"),
	 TSI = as.integer(pop == "TSI")) %>%
  select(-lab, -pop) %>%
  rename(subject = name)

write_csv(covs, "covs.csv")
