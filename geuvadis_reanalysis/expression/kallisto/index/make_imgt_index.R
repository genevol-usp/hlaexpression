devtools::load_all("/home/vitor/hlatools")
library(tidyrverse)

doMC::registerDoMC(40)

# Functions --------------------------------------------------------------------

make_index <- function(loci, infer_missing = TRUE) {

  nuc_paths <- paste0("/home/vitor/IMGTHLA/alignments/", loci, "_nuc.txt")

  if (infer_missing) {
    alignments <- 
      map_df(nuc_paths, ~hla_read_alignment(., omit = "N")) %>%
      mutate(locus = sub("^([^*]+).+$", "\\1", allele)) %>%
      select(locus, allele, cds) %>%
      group_by(locus) %>%
      filter(!all(grepl("\\*", cds))) %>%
      ungroup() %>%
      split(.$locus) %>%
      map_df(~hla_infer(., .parallel = TRUE)) %>%
      select(-locus)
  } else {
    alignments <-
      nuc_paths %>%
      map_df(~hla_read_alignment(., omit = "N", rm_incomplete = TRUE))
  }

  alignments %>%
  mutate(allele3f = ifelse(grepl("_s\\d+$", allele), allele, 
         		  hla_trimnames(allele, 3)),
         cds = hla_format_sequence(cds)) %>%
  group_by(cds) %>%
  summarize(allele = paste(allele, collapse = "/"), 
            allele3f = paste(unique(allele3f), collapse = "/")) %>%
  group_by(allele3f) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
  select(allele, cds) %>%
  arrange(allele)
}

# ------------------------------------------------------------------------------

main_loci <- c("A", "B", "C", "DQA1", "DQB1", "DRB")

other_loci <-
  list.files("/home/vitor/IMGTHLA/alignments/", pattern = "_nuc\\.txt") %>% 
  strsplit("_") %>% 
  map_chr(1) %>%
  .[! . %in% c(main_loci, "ClassI", "ClassII")]

index <- 
  bind_rows(make_index(main_loci, infer_missing = TRUE),
            make_index(other_loci, infer_missing = FALSE)) %>%
  mutate(allele = paste0("IMGT_", allele)) %>%
  split(.$allele) %>% 
  map("cds")

indexDNA <- Biostrings::DNAStringSet(unlist(index))

Biostrings::writeXStringSet(indexDNA, "./imgt_index.fa")	    
