library(GenomicAlignments)
devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- sprintf("sample_%02d", 1:50)

gencode_ids <- select(gencode_chr_tx, tx_id, gene_name)

doMC::registerDoMC(50)

bam_to_df <- function(bam) {  

    alignments <- readGAlignmentPairs(bam, use.names = TRUE) 

    tibble(name = names(alignments), 
	   reference = as.character(seqnames(alignments))) %>%
    mutate_at(vars(name, reference), function(x) sub("^([^/]+).*$", "\\1", x)) %>% 
    filter(grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", name) |
	   grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", reference)) %>%
    mutate(gene_ref = ifelse(grepl("IMGT", reference),
			     imgt_to_gname(reference),
			     tr_to_gname(reference)))
}

process_parallel <- function(paths_to_bam) {
    
    paths_to_bam %>%
	plyr::ldply(. %>% bam_to_df(), .id = "subject", .parallel = TRUE) %>%
	as_tibble() %>%
	full_join(read_ids) %>%
	mutate(read = sub("^read\\d+_", "", name),
	       gene_read = ifelse(grepl("IMGT", read), 
				  imgt_to_gname(read),
				  tr_to_gname(read))) %>%
	distinct(subject, name, reference, .keep_all = TRUE)
}

read_ids <- 
    paste0("./data/read_ids/", samples, ".txt") %>%
    setNames(samples) %>%
    plyr::ldply(. %>%
		readLines() %>%
		.[grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", .)] %>%
		sub("^@", "", .) %>%
		as_tibble, 
	    .id = "subject", .parallel = TRUE) %>%
    as_tibble() %>%
    rename(name = value) %>%
    mutate(name = sub("^([^/]+).*$", "\\1", name))

aligns_kallisto <-
    file.path("./expression/kallisto/quantifications_2", samples, "imgt.bam") %>%
    setNames(samples) %>%
    process_parallel()

aligns_star <-
    paste0("./expression/star/mappings_2/", samples, "_imgt.bam") %>%
    setNames(samples) %>%
    process_parallel()

aligns_kallisto_pri <-
    file.path("./expression/kallisto/quantifications_PRI", samples, "imgt.bam") %>%
    setNames(samples) %>%
    process_parallel()

aligns_star_pri <-
    paste0("./expression/star/mappings_PRI/", samples, "_imgt.bam") %>%
    setNames(samples) %>%
    process_parallel()

aligns_df <- 
    list(star_imgt = aligns_star, 
	 star_pri =  aligns_star_pri,
	 kallisto_imgt = aligns_kallisto, 
	 kallisto_pri = aligns_kallisto_pri) %>%
    bind_rows(.id = "aligner")

hla_genes <- paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))

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
