library(tidyverse)


link_df <- tribble(~study, ~link,
		   "geuvadis_exon", "https://doi.org/10.1038/nature12531",
		   "geuvadis_gene", "https://doi.org/10.1038/nature12531",
		   "gtex_v7", "https://doi.org/10.1038/nature24277",
		   "vince2017", "https://doi.org/10.1016/j.ajhg.2016.09.023",
		   "delaneau2018", "https://doi.org/10.1101/171694",
		   "xl9_raj2016", "https://doi.org/10.7554/eLife.12089")

read_tsv("../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/rtc/top_previous_qtls/results.tsv") %>%
    select(gene, rank, eqtl_hlapers = qtl_personalized, 
	   eqtl_previous = qtl_previous, d_prime, rtc, study) %>%
    separate_rows(study, sep = "/") %>%
    extract(study, "study", "(.+) .") %>% 
    left_join(link_df) %>%
    write_csv("./tables/S3_table.csv")
