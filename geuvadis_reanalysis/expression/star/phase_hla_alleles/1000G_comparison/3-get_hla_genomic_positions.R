devtools::load_all("~/hlaseqlib")
library(Biostrings)
library(tidyverse)

get_cds_positions <- function(nuc_file, ref_allele, cds_start, ref_cds, strand) {

    align <- hla_read_alignment(nuc_file, by_exon = TRUE) %>%
	separate(allele, c("allele", "exon"), sep = "_") %>%
	filter(allele == ref_allele) %>%
	unite(id, c("allele", "exon"), sep = "_") %>%
	mutate(cds = gsub("\\.", "", cds)) %>%
	rowwise() %>%
	mutate(cds = ifelse(strand == "+", cds, as.character(reverseComplement(DNAStringSet(cds))))) %>%
	select(id, cds) %>%
	split(.$id) %>%
	map_chr("cds") %>%
	DNAStringSet()

    align <- align[width(align) > 0]
   
    cds_pos <- matchPDict(align, DNAString(ref_cds)) %>%
	as.list() %>%
	map_df(as.data.frame, .id = 'id') %>%
	select(id, start, end) %>%
	separate(id, c("allele", "exon"), sep = "_", convert = TRUE) %>%
	as_tibble() %>%
	mutate(start = start + cds_start - 1,
	       end = end + cds_start - 1)
}

gencode_hla_exon <- 	
    get_gencode_coords("~/gencode_data/gencode.v25.annotation.gtf", "exon") %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    distinct(gene_name, strand, start, end) %>% 
    rename(locus = gene_name)

genome <- readDNAStringSet("~/gencode_data/GRCh38.primary_assembly.genome.fa.gz")
chr6 <- genome[names(genome) == "chr6 6"][[1]]
 
gen <- readDNAStringSet("/home/vitor/IMGTHLA/hla_gen.fasta")
names(gen) <- sub("^[^ ]+ ([^ ]+).*$", "\\1", names(gen))

gen_posit <- gen[grepl("^(A|DQA1|DPB1)\\*", names(gen))]
gen_negat <- reverseComplement(gen[grepl("^(B|C|DQB1|DRB1)\\*", names(gen))])

hlagen <- c(gen_posit, gen_negat)

ref_alleles <- matchPDict(hlagen, chr6) %>%
    as.list() %>%
    map_df(as.data.frame, .id = 'ref_allele') %>%
    mutate(locus = sub("^([^*]+).+$", "HLA-\\1", ref_allele)) %>%
    arrange(start) %>%
    select(locus, ref_allele, cds_start = start) %>%
    as_tibble()

call_df <- 
    tibble(ref_allele = ref_alleles$ref_allele,
	   ref_cds = as.character(hlagen[ref_alleles$ref_allele])) %>%
    left_join(ref_alleles, ., by = "ref_allele") %>%
    left_join(distinct(gencode_hla_exon, locus, strand), by ="locus") %>%
    mutate(locus_suffix = sub("HLA-", "", locus),
	   locus_suffix = ifelse(locus_suffix == "DRB1", "DRB", locus_suffix),
	   nuc_file = paste0("~/IMGTHLA/alignments/", locus_suffix, "_nuc.txt")) %>%
    select(nuc_file, ref_allele, cds_start, ref_cds, strand)

cds_positions <- pmap_df(call_df, get_cds_positions)

out <- cds_positions %>%
    mutate(locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
    left_join(gencode_hla_exon, by = "locus") %>%
    filter(start.x >= start.y, end.x <= end.y) %>%
    distinct(locus, allele, exon, start.x, end.x, .keep_all = TRUE) %>%
    select(locus, strand, allele, exon, start = start.x, end = end.x) %>%
    arrange(desc(end)) %>%
    group_by(locus) %>%
    mutate(r = min_rank(start),
	   n = 1:n()) %>%
    ungroup() %>%
    filter((strand == "+" & r >= exon) | (strand == "-" & n >= exon)) %>%
    arrange(locus, exon) %>%
    group_by(locus, allele, exon) %>%
    filter(n() == 1L) %>%
    select(locus, allele, exon, start, end)

write_tsv(out, "./cds_ref_positions.tsv")
