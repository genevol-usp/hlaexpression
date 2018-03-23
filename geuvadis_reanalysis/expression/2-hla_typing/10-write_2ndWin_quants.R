devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- readLines("../../data/sample_info/samples_phase3_ena_eur.txt")

quants <- file.path("./quantifications_top5", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    rename(allele = Name, est_counts = NumReads, tpm = TPM) %>%
    mutate(locus = imgt_to_gname(allele),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, gene_id, locus, allele, est_counts, tpm)

winners <- quants %>%
    group_by(subject, locus) %>%
    slice(which.max(est_counts)) %>%
    ungroup()

write_tsv(winners, "./quantifications_top5/winners_quants.tsv")

winners_list <- winners %>%
    select(subject, allele) %>%
    mutate(path = file.path("./quantifications_top5", subject, "winner_alleles.txt")) %>%
    split(.$subject)

map(winners_list, ~writeLines(.$allele, unique(.$path)))

non_winners_quants <- anti_join(quants, winners, by = c("subject", "allele"))

write_tsv(non_winners_quants, "./quantifications_top5/Win2_quants.tsv")
