devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name

genotypes <- 
    "../../supplemented/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("HLA-", "", locus)) %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus, allele)

write_phase_input(genotypes, "phase.inp")
