library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

sample_id <- "sample_01"
mapdir <- "mappings_2"
bam <- file.path(mapdir, sample_id, "imgt.bam")

hla_genes <- gencode_hla$gene_name
hla_regex <- "IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)"  

ids_hla_pri <- filter(gencode_pri_tx, gene_name %in% hla_genes) %>%
    pull(tx_id)

read_info <- read_tsv("../../data/read_info.tsv")

alignments <- readGAlignmentPairs(bam, use.names = TRUE) 

aligns_df <- 
    tibble(readid = names(alignments), 
	   reference = as.character(seqnames(alignments))) %>%
    mutate_at(vars(readid, reference), function(x) sub("^([^/]+).*$", "\\1", x)) %>%
    filter(readid %in% read_info$readid | grepl(hla_regex, reference)) %>%
    mutate(gene_to = ifelse(grepl("IMGT", reference),
			     imgt_to_gname(reference),
			     tr_to_gname(reference))) %>%
    full_join(read_info) %>%
    distinct(readid, reference, .keep_all = TRUE)

reads_not_aligned <- aligns_df %>%
    filter(is.na(reference)) %>%
    rowwise() %>%
    mutate(pos = list(c(m1.1:m1.2, m2.1:m2.2))) %>%
    ungroup() %>%
    select(gene_name, tx_id, pos) %>%
    unnest(pos) %>%
    count(gene_name, tx_id, pos)

write_tsv(reads_not_aligned, "./reads_not_aligned.tsv")
