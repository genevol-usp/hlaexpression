library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- tibble(subject = readLines("./data/samples.txt"),
                  code = sprintf("sample_%02d", 1:50))

gencode_ids <- select(gencode_chr_tx, tx_id, gene_name)

doMC::registerDoMC(50)

bam_to_df <- function(bam) {  

  alignments <- readGAlignmentPairs(bam, use.names = TRUE) 

  tibble(name = names(alignments), 
	 reference = as.character(seqnames(alignments))) %>%
  filter(grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", name) |
	 grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", reference)) %>%
  mutate(gene_ref = ifelse(grepl("IMGT", reference),
			   imgt_to_gname(reference),
			   tr_to_gname(reference)))
}

read_ids <- 
  paste0("./data/read_ids/", samples$code, ".txt") %>%
  setNames(samples$code) %>%
  plyr::ldply(. %>%
	      readLines() %>%
	      .[grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", .)] %>%
	      sub("^@", "", .) %>%
	      as_tibble, 
	    .id = "subject", .parallel = TRUE) %>%
  as_tibble() %>%
  rename(name = value)

aligns_kallisto <-
  file.path("./expression/kallisto/quantifications_2", samples$code, "imgt.bam") %>%
  setNames(samples$code) %>%
  plyr::ldply(. %>% bam_to_df(), .id = "subject", .parallel = TRUE) %>%
  as_tibble() %>%
  full_join(read_ids) %>%
  mutate(read = sub("^read\\d+_", "", name),
	 gene_read = ifelse(grepl("IMGT", read), 
			    imgt_to_gname(read),
			    tr_to_gname(read))) %>%
  distinct(subject, name, gene_ref, .keep_all = TRUE)

aligns_star <-
  paste0("./expression/star/mappings_2/", samples$code, "_imgt.bam") %>%
  setNames(samples$code) %>%
  plyr::ldply(. %>% bam_to_df(), .id = "subject", .parallel = TRUE) %>%
  as_tibble() %>%
  full_join(mutate(read_ids, name = sub("^([^/]+).*$", "\\1", name))) %>%
  mutate(read = sub("^read\\d+_", "", name),
	 gene_read = ifelse(grepl("IMGT", read), 
			    imgt_to_gname(read),
			    tr_to_gname(read))) %>%
  distinct(subject, name, gene_ref, .keep_all = TRUE)

aligns_df <- 
  list(kallisto = aligns_kallisto, star = aligns_star) %>%
  bind_rows(.id = "aligner")

misses <- 
  aligns_df %>%
  filter(gene_read %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  group_by(aligner, gene_read) %>%
  summarize(na = round(mean(is.na(reference)) * 100, 2)) %>%
  ungroup() %>%
  spread(aligner, na)

false_neg <-
  aligns_df %>%
  filter(gene_read %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")), !is.na(gene_ref)) %>%
  group_by(aligner, gene_read) %>%
  summarize(aver = round(mean(gene_read != gene_ref) * 100, 2)) %>%
  ungroup() %>%
  spread(aligner, aver)

false_pos <-
  aligns_df %>%
  filter(gene_ref %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")), !is.na(gene_ref)) %>%
  group_by(aligner, gene_ref) %>%
  summarize(aver = round(mean(gene_ref != gene_read) * 100, 2)) %>%
  ungroup() %>%
  spread(aligner, aver)

write_tsv(misses, "./reads_not_aligned_hla.tsv")
write_tsv(false_neg, "./reads_false_neg_hla.tsv")
write_tsv(false_pos, "./reads_false_pos_hla.tsv")
