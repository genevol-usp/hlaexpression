devtools::load_all("/home/vitor/hlatools/")
library(Biostrings)
library(tidyverse)

hla_fasta <- readDNAStringSet("./index/imgt_index.fa")
hla_fasta <- hla_fasta[grep("IMGT_(A|B|C|DQA1|DQB1|DRB1)", names(hla_fasta))] 

hla_df <-
  tibble(allele = names(hla_fasta),
	 cds = as.character(hla_fasta)) %>%
  mutate(allele = sub("IMGT_", "", allele),
	 locus = sub("^([^*]+).+$", "\\1", allele)) %>%
  select(locus, allele, cds)

ref_alleles <-
  tibble(locus = c("A", "C", "B", "DQA1", "DQB1", "DRB1"),
	 ref_allele = c("A*03:01:01", "C*07:02:01", "B*07:02:01", 
			"DQA1*01:02:01", "DQB1*06:02:01", "DRB1*15:01:01")) %>%
  left_join(hla_df, by = c("locus", "ref_allele" = "allele")) %>% 
  rename(ref_sequence = cds)

hla_df <-
  left_join(hla_df, ref_alleles, by = "locus") %>%
  select(locus, ref_allele, allele, ref_sequence, cds)

doMC::registerDoMC(6)
dist_df <-
  plyr::ddply(hla_df, ~locus, . %>% 
	      mutate(dist = map2_dbl(ref_sequence, cds, ~adist(.x, .y)/nchar(.y))),
	      .parallel = TRUE) %>%
  select(locus, allele, dist)

write_tsv(dist_df, "./distances_to_reference.tsv") 
