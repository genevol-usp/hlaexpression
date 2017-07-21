devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

genotypes <- 
  "../../quantifications_2/processed_quant.tsv" %>%
  read_tsv() %>%
  filter(locus %in% c('A', 'B', 'C', 'DQA1', 'DQB1', 'DRB1')) %>%
  mutate(subject = convert_ena_ids(subject),
         allele = gsub("IMGT_|_s\\d", "", allele),
	 allele = hla_trimnames(allele, 3)) %>%
  select(subject, locus, allele) %>%
  arrange(subject, locus, allele)

write_phase_input(genotypes, "phase.inp")
filter(genotypes, locus %in% c("A", "B", "C")) %>%
  write_phase_input("phase_class1.inp")
filter(genotypes, locus %in% c("DQA1", "DQB1", "DRB1")) %>%
  write_phase_input("phase_class2.inp")
