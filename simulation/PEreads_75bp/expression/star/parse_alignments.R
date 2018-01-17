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

read_summary <- aligns_df %>%
    filter(gene_from %in% hla_genes) %>%
    group_by(name) %>%
    summarize(gene_from = unique(gene_from),
	      not_aligned = as.integer(all(is.na(gene_to))),
	      maps_to_original = as.integer(any(!is.na(gene_to) & gene_to == gene_from)),
	      n_genes_to = n_distinct(gene_to),
	      n_refs_to = n_distinct(reference)) %>%
    group_by(gene_from) %>%
    summarize(perc_not_aligned = mean(not_aligned) * 100,
	      perc_not_aligned_to_original = mean(not_aligned == 0 & maps_to_original == 0) * 100,
	      perc_aligned_to_original_multimap = mean(not_aligned == 0 & maps_to_original == 1 & n_genes_to > 1) * 100,
	      perc_aligned_to_original_uniquely = mean(not_aligned == 0 & maps_to_original == 1 & n_genes_to == 1) * 100) %>%
    ungroup()

alignments_to_diff_gene <- aligns_df %>%
    filter(gene_from %in% hla_genes, !is.na(gene_to)) %>%
    count(gene_from, gene_to) %>% 
    group_by(gene_from) %>% 
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    filter(gene_from != gene_to) %>%
    select(gene_from, gene_to, perc)
    
alignments_from_diff_gene <- aligns_df %>%
    filter(gene_to %in% hla_genes, !is.na(gene_to)) %>%
    count(gene_from, gene_to) %>%
    group_by(gene_to) %>%
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    filter(gene_from != gene_to) %>%
    select(gene_to, gene_from, perc)

write_tsv(read_summary, 
	  file.path(mapdir, sample_id, "read_summary_hla.tsv"))

write_tsv(alignments_to_diff_gene, 
	  file.path(mapdir, sample_id, "alignments_to_diff_gene_hla.tsv"))

write_tsv(alignments_from_diff_gene, 
	  file.path(mapdir, sample_id, "alignments_from_diff_gene_hla.tsv"))
