devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

genotypes <- 
    "../../quantifications_2/processed_quant.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("HLA-", "", locus)) %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus, allele)

write_phase_input(genotypes, "phase.inp")
