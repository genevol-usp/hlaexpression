devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

samples <- geuvadis_info %>% 
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications/processed_imgt_quants.tsv")

calls <- imgt_quants %>%
    filter(locus %in% hla_genes) %>%
    select(subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(as.character(subject)),
	   locus = sub("^HLA-", "", locus),
	   allele = hla_trimnames(sub("IMGT_", "", allele), 3)) %>%
    arrange(subject, locus, allele)

typing_errors <- 
    calc_genotyping_accuracy(calls, pag, by_locus = FALSE) %>%
    group_by(subject, locus) %>%
    filter(any(!correct)) %>%
    ungroup()

manually_accepted <-
    tribble(~subject, ~locus,
	    "HG00115", "A",
	    "HG00150", "DQB1",
	    "HG00246", "A",
	    "HG00256", "DQB1",
	    "HG00273", "A",
	    "HG00311", "A",
	    "HG00360", "B",
	    "NA06986", "DQB1",
	    "NA10847", "C",
	    "NA11830", "DQB1",
	    "NA11832", "B",
	    "NA11840", "B",
	    "NA11892", "DQB1",
	    "NA11920", "DRB1",
	    "NA11930", "DQB1",
	    "NA11994", "A",
	    "NA12005", "B",
	    "NA12045", "DRB1",
	    "NA12058", "B",
	    "NA12058", "C",
	    "NA12156", "DQB1",
	    "NA12234", "DQB1",
	    "NA12272", "DRB1",
	    "NA12286", "DQB1",
	    "NA20534", "A")

typing_errors_out <- anti_join(typing_errors, manually_accepted) %>%
    distinct(subject, locus)

write_tsv(typing_errors_out, "./typing_errors.tsv")
