devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(Biostrings)
library(tidyverse)

hla_genes <- c(gencode_hla$gene_name, "HLA-E", "HLA-G")

gencode_hla2 <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes)

hla_positions <- read_tsv("./cds_ref_positions.tsv") %>%
    select(-allele)

hla_cds <- read_tsv("./index_ref_positions.tsv") %>%
	mutate(locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
	inner_join(hla_positions, by = c("locus", "exon")) %>%
	left_join(select(gencode_hla2, locus = gene_name, strand), by = "locus")

hla_cds_neg <- hla_cds %>%
    filter(strand == "-") %>%
    rowwise() %>%
    mutate(cds = as.character(reverseComplement(DNAString(cds)))) %>%
    ungroup()

hla_cds_final <- hla_cds %>%
    filter(strand == "+") %>%
    bind_rows(hla_cds_neg) %>%
    select(-locus, -strand) %>%
    arrange(allele, exon)

hla_cds_pos_df <- hla_cds_final %>%
    mutate(cds = strsplit(cds, ""), 
	   pos = map2(start, end, `:`)) %>%
    select(-start, -end) %>%
    unnest()

genos <- 
    "../../expression/2-hla_typing/genotype_calls.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    left_join(geuvadis_info, by = c("subject" = "ena_id")) %>%
    select(subject = name, locus, allele) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    arrange(subject, locus, allele)

# Some alleles have the same CDS when only non-del positions in the
# HLA reference allele are considered, but have different CDS when the whole
# sequence is considered. For example, DQB1*05:03:01:01 and DQB1*05:03:01:03
genos_fix <- genos %>%
    mutate(allele_3F = hla_trimnames(allele, 3)) %>%
    left_join(distinct(hla_cds, locus, allele), by = "allele") %>%
    mutate(allele_fix = ifelse(is.na(locus.y), allele_3F, allele)) %>%
    select(subject, locus = locus.x, allele = allele_fix)

hla_genos_and_sequence <- left_join(genos_fix, hla_cds_pos_df, by = "allele")

vcf <- read_tsv("./hla_snps.vcf", comment = "##", progress = FALSE) %>%
    filter(nchar(REF) == 1L, grepl("^[ACGT](,[ACGT])*$", ALT)) %>%
    select(-1, -3, -(6:9)) %>%
    gather(subject, genotype, -(1:3)) %>%
    filter(subject %in% genos$subject) %>%
    separate(genotype, c("h1", "h2"), sep = "\\|", convert = TRUE, remove = FALSE) %>%
    mutate(ALT = gsub(",", "", ALT),
	   h1 = ifelse(h1 == 0, REF, substring(ALT, h1, h1)),
	   h2 = ifelse(h2 == 0, REF, substring(ALT, h2, h2))) %>%
    select(subject, pos = POS, genotype, h1, h2)

out_df <-
    inner_join(vcf, hla_genos_and_sequence, by = c("subject", "pos")) %>%
    group_by(subject, locus, pos) %>%
    filter(!any(cds == "N")) %>% #for a given pair of alleles, consider only sequenced positions
    ungroup() %>%
    mutate(a1 = cds == h1, 
	   a2 = cds == h2) %>%
    group_by(subject, locus, exon, allele) %>%
    summarise(a1 = sum(!a1), a2 = sum(!a2)) %>%
    ungroup() %>%
    gather(hap, diffs, a1, a2) %>%
    mutate(hap = as.integer(sub("^a", "", hap))) %>%
    arrange(subject, locus, exon, allele, hap)

out_df %>%
    group_by(subject, locus, allele, hap) %>%
    summarise(diffs = sum(diffs)) %>%
    write_tsv("./hla_diffs_to_1000Ghaps.tsv")

write_tsv(out_df, "./hla_diffs_to_1000Ghaps_by_exon.tsv")
