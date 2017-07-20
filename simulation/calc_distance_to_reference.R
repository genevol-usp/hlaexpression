devtools::load_all("/home/vitor/hlaseqlib")
library(Biostrings)
library(tidyverse)

hla_fasta <- readDNAStringSet("../geuvadis_reanalysis/expression/kallisto/index/imgt_index.fa")
hla_fasta <- hla_fasta[grep("IMGT_(A|B|C|DQA1|DQB1|DRB1)", names(hla_fasta))] 

hla_df <-
  tibble(allele = names(hla_fasta),
	 cds = as.character(hla_fasta)) %>%
  mutate(allele = sub("IMGT_", "", allele),
	 locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
  select(locus, allele, cds)

ref_alleles <- 
  "../geuvadis_reanalysis/expression/kallisto/phase_hla_alleles/1000G_comparison/hla_ref_alleles.tsv" %>%
  read_tsv() %>%
  mutate(allele = hla_trimnames(allele, 3)) %>%
  left_join(hla_df, by = c("locus", "allele")) %>% 
  rename(ref_allele = allele, ref_sequence = cds)

hla_df <-
  left_join(hla_df, ref_alleles, by = "locus") %>%
  select(locus, ref_allele, allele, ref_sequence, cds)

dist_d <-
  hla_df %>%
  split(.$locus) %>%
  parallel::mclapply(. %>% 
    mutate(dist = map2_dbl(ref_sequence, cds, ~adist(.x, .y)/nchar(.y))),
    mc.cores = 6) %>%
  bind_rows() %>%
  select(locus, allele, dist)

write_tsv(dist_df, "./distances_to_reference.tsv") 
