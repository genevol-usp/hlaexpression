devtools::load_all("/home/vitor/hlaseqlib")
library(Biostrings)
library(tidyverse)

map_to_genome <- function(locus, strand, annot_start, annot_end, ref_allele, chr6) {

    locus <- sub("HLA-", "", locus) 
    if (locus == "DRB1") locus <- "DRB"
   
    ref_seq <- 
	paste0("/home/vitor/IMGTHLA/alignments/", locus, "_nuc.txt") %>%
	hla_read_alignment(omit = "N", rm_incomplete = TRUE, by_exon = TRUE) %>%
	separate(allele, c("allele", "exon"), sep = "_") %>%
	filter(allele == ref_allele) %>%
	mutate(cds = hla_format_sequence(cds)) %>%
	filter(cds != "") %>%
	select(allele, exon, cds) %>%
	unite(allele, allele:exon) %>%
	group_by(allele) %>%
	mutate(cds = ifelse(strand == "+", cds, as.character(reverseComplement(DNAStringSet(cds))))) %>%
	ungroup() %>%
	select(allele, cds)

    x_seqs <- 
	ref_seq %>%
	split(.$allele) %>%
	map_chr("cds") %>%
	DNAStringSet()

    matchPDict(x_seqs, chr6) %>%
	as.list() %>%
	map_df(. %>% as.data.frame, .id = "allele_exon") %>%
	filter(start %in% unlist(annot_start) | end %in% unlist(annot_end)) %>%
	arrange(start) %>%
	group_by(allele_exon) %>% 
	mutate(coords = list(start:end)) %>%
	ungroup() %>%
	summarise(coords = list(unlist(coords))) %>%
	pull(coords)
}

genome <- 
    "/home/vitor/gencode_data/GRCh38.primary_assembly.genome.fa.gz" %>%
    readDNAStringSet()

chr6 <- genome[names(genome) == "chr6 6"][[1]]

samples <- 
    read_tsv("../../../../data/sample_info/samples_phase3.tsv") %>%
    pull(subject)

ref_alleles <- read_tsv("../../../../../imgt_index/hla_ref_alleles.tsv")

genos <- 
    read_tsv("../../quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
    select(subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(subject),
	   allele = sub("^([^=]+).*$", "\\1", allele))

vcf <-
    read_tsv("./hla_snps.vcf", comment = "##", progress = FALSE) %>%
    filter(nchar(REF) == 1, grepl("^[ACGT](,[ACGT])*$", ALT)) %>% 
    select(-1, -3, -6, -7, -8, -9) %>%
    gather(subject, genotype, HG00096:NA21144) %>%
    filter(subject %in% samples) %>%
    separate(genotype, c("h1", "h2"), sep = "\\|") %>%
    mutate(ALT = gsub(",", "", ALT),
	   h1 = ifelse(h1 == 0, REF, substring(ALT, h1, h1)),
	   h2 = ifelse(h2 == 0, REF, substring(ALT, h2, h2))) %>%
    select(subject, pos = POS, h1, h2)

duplicated_pos <- vcf %>%
    count(subject, pos) %>%
    distinct(pos, n) %>%
    filter(n > 1L) %>%
    pull(pos)

vcf <- vcf %>% filter(! pos %in% duplicated_pos) 

hla_db <-
    "/home/vitor/gencode_data/gencode.v25.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "exon") %>%
    select(locus = gene_name, strand, start, end) %>%
    filter(locus %in% ref_alleles$locus) %>%
    distinct() %>%
    group_by(locus, strand) %>%
    summarize(annot_start = list(start), annot_end = list(end)) %>%
    left_join(ref_alleles, by = "locus") %>%
    ungroup() %>%
    rename(ref_allele = allele)

pos_df <- 
    hla_db %>%
    mutate(pos = pmap(hla_db, map_to_genome, chr6 = chr6),
	   pos = flatten(pos)) %>%
    select(locus, strand, pos)
  
allele_index <- 
    readDNAStringSet("../../../../../imgt_index/index_ref_positions.fa")

allele_seq_df <- 
    tibble(allele = names(allele_index), cds = as.character(allele_index)) %>%
    mutate(locus = imgt_to_gname(allele)) %>%
    select(locus, allele, cds) %>%
    filter(allele %in% genos$allele) %>%
    left_join(pos_df, by = "locus") %>%
    group_by(allele) %>%
    mutate(cds = map_chr(cds, ~ifelse(strand == "+", ., as.character(reverseComplement(DNAStringSet(.)))))) %>%
    ungroup() %>%
    mutate(cds = strsplit(cds, "")) %>%
    unnest(pos, cds) %>%
    filter(pos %in% vcf$pos) %>%
    select(locus, allele, pos, snp_allele = cds)

variant_df <-
    left_join(genos, allele_seq_df, by = c("locus", "allele")) %>%
    left_join(vcf, by = c("subject", "pos")) %>%
    mutate(a1 = snp_allele == h1,	a2 = snp_allele == h2) %>%
    group_by(subject, locus, allele) %>%
    summarize(a1 = sum(!a1), a2 = sum(!a2)) %>%
    ungroup() %>%
    gather(hap, diffs, a1, a2) %>%
    arrange(subject, allele, hap)

write_tsv(variant_df, "./hla_haps_mapped_to_1000G.tsv")
