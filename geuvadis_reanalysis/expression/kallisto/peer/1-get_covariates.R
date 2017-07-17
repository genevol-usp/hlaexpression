devtools::load_all("/home/vitor/genomicRutils")

samples <- readr::read_lines("../samples_phase3.txt")

labs <- levels(reorder(geuvadis_info$lab, geuvadis_info$lab_code))
pops <- unique(geuvadis_info$pop)

matrix_lab <- 
  matrix(0, nrow = length(samples), ncol = length(labs), 
	 dimnames = list(NULL, labs)) %>%
  tibble::as_tibble()

matrix_pop <- 
  matrix(0, nrow = length(samples), ncol = length(pops), 
	 dimnames = list(NULL, pops)) %>%
  tibble::as_tibble()

covs <- 
  geuvadis_info %>%
  dplyr::filter(name %in% samples) %>%
  dplyr::select(name, lab, pop) %>%
  dplyr::bind_cols(matrix_lab, matrix_pop) %>%
  dplyr::mutate(UNIGE = as.integer(lab == "UNIGE"),
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
  dplyr::select(-lab, -pop) %>%
  dplyr::rename(subject = name)

readr::write_csv(covs, "covs.csv")
