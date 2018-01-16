library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

sample_id <- commandArgs(TRUE)[1]
mapdir <- commandArgs(TRUE)[2]
bam <- file.path(mapdir, sample_id, "imgt.bam")

hla_genes <- gencode_hla$gene_name
hla_regex <- "IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)"  

ids_hla_pri <- readLines("../../../../imgt_index/hla_ids_pri.txt")

read_ids <- paste0("../../data/read_ids/", sample_id, ".txt") %>%
    readLines() %>%
    .[grepl(hla_regex, .)] %>%
    sub("^@", "", .) %>%
    as_tibble() %>%
    rename(name = value) %>%
    mutate(name = sub("^([^/]+).*$", "\\1", name))

alignments <- readGAlignmentPairs(bam, use.names = TRUE) 

aligns_df <- 
    tibble(name = names(alignments), 
	   reference = as.character(seqnames(alignments))) %>%
    mutate_at(vars(name, reference), function(x) sub("^([^/]+).*$", "\\1", x)) %>% 
    filter((name %in% ids_hla_pri | grepl("IMGT", name)) | 
	   (reference %in% ids_hla_pri | grepl("IMGT", reference))) %>% 
    mutate(gene_to = ifelse(grepl("IMGT", reference),
			     imgt_to_gname(reference),
			     tr_to_gname(reference))) %>%
    full_join(read_ids) %>%
    mutate(read = sub("^read\\d+_", "", name),
	   gene_from = ifelse(grepl("IMGT", read), 
			      imgt_to_gname(read),
			      tr_to_gname(read))) %>%
    distinct(name, reference, .keep_all = TRUE)

reads_not_aligned <- aligns_df %>%
    filter(gene_from %in% hla_genes) %>%
    group_by(gene_from) %>%
    summarize(perc = round(mean(is.na(reference)) * 100, 2)) %>%
    ungroup()

reads_lost_to_other_genes <- aligns_df %>%
    filter(gene_from %in% hla_genes, !is.na(gene_to)) %>%
    count(gene_from, gene_to) %>% 
    group_by(gene_from) %>% 
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    filter(gene_from != gene_to) %>%
    select(gene_from, gene_to, perc)
    
reads_gained_from_other_genes <- aligns_df %>%
    filter(gene_to %in% hla_genes, !is.na(gene_to)) %>%
    count(gene_from, gene_to) %>%
    group_by(gene_to) %>%
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    filter(gene_from != gene_to) %>%
    select(gene_to, gene_from, perc)

write_tsv(reads_not_aligned, 
	  file.path(mapdir, sample_id, "reads_not_aligned_hla.tsv"))

write_tsv(reads_lost_to_other_genes, 
	  file.path(mapdir, sample_id, "reads_lost_to_other_genes_hla.tsv"))

write_tsv(reads_gained_from_other_genes, 
	  file.path(mapdir, sample_id, "reads_gained_from_other_genes_hla.tsv"))
