devtools::load_all("/home/vitor/hlaseqlib")
library(Biostrings)
library(tidyverse)

main_loci <- c("A", "B", "C", "DQA1", "DQB1", "DRB")

other_loci <-
  list.files("/home/vitor/IMGTHLA/alignments/", pattern = "_nuc\\.txt") %>% 
  strsplit("_") %>% 
  map_chr(1) %>%
  .[! . %in% c(main_loci, "ClassI", "ClassII")]

loci_df <- 
  tibble(locus = c(main_loci, other_loci)) %>%
  mutate(infer = locus %in% main_loci,
	 seqs = map2(locus, infer, hla_make_sequences, n_cores = 16)) 

seqs_df <- select(loci_df, seqs) %>% unnest()

index_df <- 
  seqs_df %>%
  mutate(cds = hla_format_sequence(cds)) %>%
  group_by(cds) %>%
  summarize(allele = paste(allele, collapse = "/")) %>%
  ungroup() %>%
  mutate(allele3f = hla_trimnames(allele)) %>%
  group_by(allele3f) %>%
  mutate(n = n()) %>%
  ungroup() %>% 
  arrange(allele)

index <-
  index_df %>%
  mutate(allele = paste0("IMGT_", ifelse(n > 1L, allele, allele3f))) %>%
  select(allele, cds) %>%
  split(.$allele) %>%
  map_chr("cds") %>%
  DNAStringSet()

writeXStringSet(index, "./imgt_index.fa")	    

index_info <-
  index_df %>%
  select(allele, allele3f, n) %>%
  rownames_to_column() %>%
  separate_rows(allele, sep = "/")

ref_pos_df <- 
  read_tsv("../phase_hla_alleles/1000G_comparison/hla_ref_alleles.tsv") %>%
  left_join(seqs_df, by = "allele") %>%
  mutate(pos = map(cds, ~which(unlist(strsplit(., "")) != "."))) %>%
  select(locus, pos)

index_ref_pos_df <-
  seqs_df %>%
  filter(grepl("^(A|B|C|DQA1|DQB1|DRB1)", allele)) %>%
  mutate(locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
  left_join(ref_pos_df, by = "locus") %>%
  mutate(cds = map2_chr(cds, pos, ~paste(substring(.x, .y, .y), collapse = "")),
	 cds = gsub("\\.", "-", cds),
	 cds = gsub("\\*", "N", cds)) %>%
  select(allele, cds) %>%
  left_join(index_info, by = "allele") %>%
  group_by(rowname) %>%
  summarize(allele = paste(allele, collapse = "/"),
	    allele3f = unique(allele3f),
	    n = unique(n),
	    cds = unique(cds)) %>%
  ungroup() %>%
  arrange(allele)

index_ref_pos <-
  index_ref_pos_df %>%
  mutate(allele = paste0("IMGT_", ifelse(n > 1L, allele, allele3f))) %>%
  select(allele, cds) %>%
  split(.$allele) %>%
  map_chr("cds") %>%
  DNAStringSet()

writeXStringSet(index_ref_pos, "./index_ref_positions.fa")
