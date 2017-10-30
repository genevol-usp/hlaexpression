library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

sample_id <- commandArgs(TRUE)[1]
bam <- commandArgs(TRUE)[2]
outdir <- commandArgs(TRUE)[3]

gencode_ids <- select(gencode_chr_tx, tx_id, gene_name)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))
hla_regex <- "IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)"  

ids_hla_pri <- readLines("../../../../imgt_index/hla_ids_pri.txt")

read_ids <- 
    paste0("../../data/read_ids/", sample_id, ".txt") %>%
    readLines() %>%
    .[grepl(hla_regex, .)] %>%
    sub("^@", "", .) %>%
    as_tibble() %>%
    rename(name = value) %>%
    mutate(subject = sample_id,
	   name = sub("^([^/]+).*$", "\\1", name))

alignments <- readGAlignmentPairs(bam, use.names = TRUE) 

aligns_df <- 
    tibble(name = names(alignments), 
	   reference = as.character(seqnames(alignments))) %>%
    mutate_at(vars(name, reference), function(x) sub("^([^/]+).*$", "\\1", x)) %>% 
    filter((name %in% ids_hla_pri | grepl("IMGT", name)) | 
	   (reference %in% ids_hla_pri | grepl("IMGT", reference))) %>% 
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
    inner_join(aligns_df, read_ids) %>%
    group_by(gene_read) %>%
    summarize(na = round(mean(is.na(reference)) * 100, 2)) %>%
    ungroup()

aligned_to_diff <-
    inner_join(aligns_df, read_ids) %>%
    filter(!is.na(gene_ref)) %>%
    group_by(gene_read) %>%
    summarize(perc = round(mean(gene_read != gene_ref) * 100, 2)) %>%
    ungroup()

aligned_summary <-
    inner_join(aligns_df, read_ids) %>%
    filter(!is.na(gene_ref)) %>%
    group_by(gene_read, gene_ref) %>%
    summarize(n = n()) %>%
    group_by(gene_read) %>%
    mutate(perc = round(n/sum(n) * 100, 4)) %>%
    arrange(gene_read, desc(perc))

false_pos <-
    aligns_df %>%
    filter(gene_ref %in% hla_genes, !is.na(gene_ref)) %>%
    group_by(gene_ref) %>%
    summarize(perc = round(mean(gene_ref != gene_read) * 100, 2)) %>%
    ungroup()

write_tsv(not_aligned, paste0(outdir, "/reads_not_aligned_hla.tsv"))
write_tsv(aligned_to_diff, paste0(outdir, "/reads_false_neg_hla.tsv"))
write_tsv(aligned_summary, paste0(outdir, "/hla_aligned_summary.tsv"))
write_tsv(false_pos, paste0(outdir, "/reads_false_pos_hla.tsv"))
