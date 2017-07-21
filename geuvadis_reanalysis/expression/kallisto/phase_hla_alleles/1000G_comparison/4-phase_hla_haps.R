devtools::load_all("~/genomicRutils")

doMC::registerDoMC(48)

hla_haps <- readr::read_tsv("./hla_haps_mapped_to_1000G.tsv")

phased <- plyr::ddply(hla_haps, ~subject, . %>% phase_hla, 
		      .id = "subject", .parallel = TRUE)

readr::write_tsv(phased, "./hla_haps_mapped_to_1000G_phased.tsv")
