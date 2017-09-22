library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

bam_path <- commandArgs(TRUE)[1]

gencode_ids <- select(gencode_chr_tx, tx_id, gene_name)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))
hla_regex <- "IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)"  

ids_to_filter <- readLines("../data/ids_to_filter.txt")

read_ids <- 
    paste0("../data/read_ids/", sample_id, ".txt") %>%
    readLines() %>%
    .[grepl(hla_regex, .)] %>%
    sub("^@", "", .) %>%
    as_tibble() %>%
    rename(name = value) %>%
    mutate(subject = sample_id,
	   name = sub("^([^/]+).*$", "\\1", name))

alignments <- readGAlignmentPairs(bam_path, use.names = TRUE) 

aligns_df <- 
    tibble(name = names(alignments), 
	   reference = as.character(seqnames(alignments))) %>%
    mutate_at(vars(name, reference), function(x) sub("^([^/]+).*$", "\\1", x)) %>% 
    filter(name %in% ids_to_filter | reference %in% ids_to_filter) %>% 
    mutate(gene_ref = ifelse(grepl("IMGT", reference),
			     imgt_to_gname(reference),
			     tr_to_gname(reference))) %>%
    full_join(read_ids) %>%
    mutate(read = sub("^read\\d+_", "", name),
	   gene_read = ifelse(grepl("IMGT", read), 
			      imgt_to_gname(read),
			      tr_to_gname(read))) %>%
    distinct(subject, name, reference, .keep_all = TRUE)

not_aligned <- 
    aligns_df %>%
    filter(gene_read %in% hla_genes) %>%
    group_by(aligner, gene_read) %>%
    summarize(na = round(mean(is.na(reference)) * 100, 2)) %>%
    ungroup() %>%
    spread(aligner, na)

aligned_to_diff <-
    aligns_df %>%
    filter(gene_read %in% hla_genes, !is.na(gene_ref)) %>%
    group_by(aligner, gene_read) %>%
    summarize(aver = round(mean(gene_read != gene_ref) * 100, 2)) %>%
    ungroup() %>%
    spread(aligner, aver)

diff_refs_df <-
    aligns_df %>%
    filter(gene_read %in% hla_genes, !is.na(gene_ref), 
	   gene_read != gene_ref) %>%
    count(aligner, gene_read, gene_ref)

false_pos <-
    aligns_df %>%
    filter(gene_ref %in% hla_genes, !is.na(gene_ref)) %>%
    group_by(aligner, gene_ref) %>%
    summarize(aver = round(mean(gene_ref != gene_read) * 100, 2)) %>%
    ungroup() %>%
    spread(aligner, aver)

false_refs_df <-
    aligns_df %>%
    filter(gene_read %in% hla_genes, !is.na(gene_ref), 
	   gene_ref != gene_read) %>%
    count(aligner, gene_read, gene_ref)

write_tsv(not_aligned, "./reads_not_aligned_hla.tsv")
write_tsv(aligned_to_diff, "./reads_false_neg_hla.tsv")
write_tsv(false_pos, "./reads_false_pos_hla.tsv")
write_tsv(diff_refs_df, "./diff_refs_hla.tsv")
write_tsv(false_refs_df, "./false_refs_hla.tsv")
