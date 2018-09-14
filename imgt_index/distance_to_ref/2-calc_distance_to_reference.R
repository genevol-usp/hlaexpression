devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(Biostrings)
library(tidyverse)

index <- readDNAStringSet("../imgt_index.fa")

hla_regex <- sub("HLA-", "", gencode_hla$gene_name) %>%
    paste(collapse = "|") %>%
    paste0("IMGT_(", ., ")")

hla_fasta <- index[grep(hla_regex, names(index))] 

hla_df <-
    tibble(allele = names(hla_fasta), cds = as.character(hla_fasta)) %>%
    mutate(allele = sub("IMGT_", "", allele),
	   locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    select(locus, allele, cds)

ref_alleles <- read_tsv("./ref_alleles.tsv") %>%
    mutate(ref_allele = hla_trimnames(ref_allele, 3)) %>%
    left_join(hla_df, by = c("locus", "ref_allele" = "allele")) %>% 
    rename(ref_sequence = cds)

hla_df <-
    left_join(hla_df, ref_alleles, by = "locus") %>%
    select(locus, ref_allele, allele, ref_sequence, cds)

dist_df <- hla_df %>%
    mutate(dist = map2_dbl(cds, ref_sequence, ~adist(.x, .y)/nchar(.y))) %>%
    select(locus, allele, dist)

write_tsv(dist_df, "./distances_to_reference.tsv") 
