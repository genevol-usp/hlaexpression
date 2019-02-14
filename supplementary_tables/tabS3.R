library(tidyverse)

link_df <- tribble(~study, ~link,
		   "GeuvadisExon", "https://doi.org/10.1038/nature12531",
		   "GeuvadisGene", "https://doi.org/10.1038/nature12531",
		   "GTExV7", "https://doi.org/10.1038/nature24277",
		   "Fairfax2012", "https://doi.org/10.1038/ng.2205")

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/top_previous_qtls/results.tsv") %>%
    select(gene, rank, eqtl_hlapers = qtl_personalized, 
	   eqtl_previous = qtl_previous, d_prime, rtc, info) %>%
    extract(info, "study", "(.+) .") %>% 
    left_join(link_df) %>%
    write_csv("./tables/S3_table.csv")
