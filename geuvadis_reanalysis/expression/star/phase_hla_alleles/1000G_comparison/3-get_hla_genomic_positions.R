library(Biostrings)
library(tidyverse)

genome <- readDNAStringSet("~/gencode_data/GRCh38.primary_assembly.genome.fa.gz")
chr6 <- genome[names(genome) == "chr6 6"][[1]]
 
gen <- readDNAStringSet("/home/vitor/IMGTHLA/hla_gen.fasta")
gen_posit <- gen[grepl("^[^ ]+ (A|DQA1|DPB1)\\*", names(gen))]
gen_negat <- reverseComplement(gen[grepl("^[^ ]+ (B|C|DQB1|DRB1)\\*", names(gen))])

hlagen <- c(gen_posit, gen_negat)

ref_alleles <- 
  matchPDict(hlagen, chr6) %>%
  as.list() %>%
  map_df(as.data.frame, .id = 'allele') %>%
  mutate(allele = map_chr(strsplit(allele, " "), 2),
	 locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
  arrange(start) %>%
  select(locus, allele)

write_tsv(ref_alleles, "./hla_ref_alleles.tsv")
