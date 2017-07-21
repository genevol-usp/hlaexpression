devtools::load_all("/home/vitor/hlaseqlib")
library(Biostrings)
library(tidyverse)

map_to_genome <- function(locus, ref_allele, strand, annot_starts, annot_ends, chr6seq) {

  locus <- sub("HLA-", "", locus) 
  if (locus == "DRB1") locus <- "DRB"
  
  nuc <- paste0("/home/vitor/IMGTHLA/alignments/", locus, "_nuc.txt")

  ref_seq <- 
    hla_read_alignment(nuc, omit = "N", rm_incomplete = TRUE, by_exon = TRUE) %>%
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
  
  x_seqs <- DNAStringSet(ref_seq$cds)
  names(x_seqs) <- ref_seq$allele

  genomicpos <- 
    matchPDict(x_seqs, chr6seq) %>%
    as.list() %>%
    map_df(. %>% as.data.frame, .id = "allele_exon") %>%
    filter(start %in% unlist(annot_starts) | end %in% unlist(annot_ends)) %>%
    arrange(start) %>%
    group_by(allele_exon) %>% 
    mutate(coords = list(start:end)) %>%
    ungroup() %>%
    summarise(coords = list(unlist(coords))) %>%
    pull(coords)

  alignments <- hla_read_alignment(nuc, omit = "N")
  
  if (locus == "DRB") {
    alignments <- filter(alignments, !grepl("^DRB[2-9]", allele))
  }

  locus_pos <-
    which(unlist(strsplit(alignments$cds[alignments$allele == ref_allele], "")) != ".")
   
  alignments %>%
  mutate(allele = hla_trimnames(allele, 3),
	 pos = genomicpos, 
	 cds = map_chr(cds, ~paste(unlist(strsplit(., ""))[locus_pos], collapse = "")),
	 cds = gsub("\\.", "-", cds),
	 cds = gsub("\\*", "N", cds),
	 cds = map_chr(cds, ~ifelse(strand == "+", ., as.character(reverseComplement(DNAStringSet(.)))))) %>%
  distinct(allele, cds, .keep_all = TRUE) %>%
  mutate(cds = strsplit(cds, "")) %>%
  unnest(pos, cds) %>%
  rename(snp_allele = cds)
}

genome <- readDNAStringSet("~/gencode_data/GRCh38.primary_assembly.genome.fa.gz")
chr6 <- genome[names(genome) == "chr6 6"][[1]]

samples <- 
  read_tsv("../../samples_phase3.tsv") %>%
  pull(subject)

ref_alleles <- read_tsv("./hla_ref_alleles.tsv")

genos <- 
  read_tsv("../../quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1")) %>%
  select(subject, locus, allele) %>%
  mutate(subject = convert_ena_ids(subject),
	 allele = gsub("IMGT_|_s\\d", "", allele),
	 allele = hla_trimnames(allele, 3))

vcf <-
  read_tsv("./hla_snps.vcf", comment = "##") %>%
  filter(nchar(REF) == 1, grepl("^[ACGT](,[ACGT])*$", ALT)) %>%
  select(-1, -3, -6, -7, -8, -9) %>%
  gather(subject, genotype, HG00096:NA21144) %>%
  filter(subject %in% samples) %>%
  separate(genotype, c("h1", "h2"), sep = "\\|") %>%
  mutate(ALT = gsub(",", "", ALT),
         h1 = ifelse(h1 == 0, REF, substring(ALT, h1, h1)),
         h2 = ifelse(h2 == 0, REF, substring(ALT, h2, h2))) %>%
  select(subject, pos = POS, h1, h2)

hla_db <-
  get_gencode_coords("~/gencode_data/gencode.v26.annotation.gtf.gz", feature = "exon") %>%
  select(locus = gene_name, strand, start, end) %>%
  filter(locus %in% ref_alleles$locus) %>%
  distinct() %>%
  group_by(locus, strand) %>%
  summarize(start = list(start), end = list(end)) %>%
  left_join(ref_alleles, by = "locus") 

out <-
  map_df(split(hla_db, hla_db$locus), 
	 ~map_to_genome(.x$locus, .x$allele, .x$strand, .x$start, .x$end, chr6)) %>% 
  filter(pos %in% vcf$pos)

locus_pos <- out %>%
  mutate(locus = sub("^([^*]+).+$", "\\1", allele)) %>%
  distinct(locus, pos)

vcf <- 
  inner_join(vcf, locus_pos, by = "pos") %>%
  inner_join(genos, by = c("subject", "locus"))

variant_df <- 
  inner_join(vcf, out, by = c("allele", "pos")) %>%
  mutate(a1 = snp_allele == h1 | snp_allele == "N",
	 a2 = snp_allele == h2 | snp_allele == "N") %>%
  group_by(subject, locus, allele) %>%
  summarize(a1 = sum(!a1), a2 = sum(!a2)) %>%
  ungroup() %>%
  gather(hap, diffs, a1, a2) %>%
  arrange(subject, hap, allele)

write_tsv(variant_df, "./hla_haps_mapped_to_1000G.tsv")
