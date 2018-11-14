devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

locus <- commandArgs(TRUE)[1]

if (grepl("DRB\\d+", locus)) { 
    locus_nuc <- "DRB"
} else {
    locus_nuc <- locus
}

nuc_file <- paste0("/home/vitor/IMGTHLA/alignments/", locus_nuc, "_nuc.txt")

hla_df <- hla_read_alignment(nuc_file, by_exon = TRUE) %>%
    mutate(gene = sub("^([^*]+).+$", "\\1", allele)) %>%
    filter(gene == locus) %>%
    select(-gene) %>%
    separate(allele, c("allele", "exon"), sep = "_") %>%
    group_by(allele) %>%
    summarise(cds = paste(cds, collapse = "|")) %>%
    ungroup()

ref_seq <- read_tsv("./cds_ref_positions.tsv") %>%
    distinct(locus, allele) %>%
    rename(gene = locus) %>%
    filter(gene == paste0("HLA-", locus)) %>%
    left_join(hla_df, by = "allele") %>%
    pull(cds) %>%
    strsplit("") %>%
    unlist() 

ref_pos <- which(! ref_seq == ".")

hla_df$cds <- str_split(hla_df$cds, "", simplify = TRUE) %>%
    .[, ref_pos] %>%
    apply(1, . %>% paste(collapse = ""))

out_df <- hla_df %>%
    mutate(allele3f = hla_trimnames(allele, 3)) %>%
    distinct(allele3f, cds, .keep_all = TRUE) %>%
    group_by(allele3f) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
    mutate(cds = strsplit(cds, "\\|")) %>%
    unnest() %>%
    group_by(allele) %>%
    mutate(exon = seq_len(n())) %>%
    filter(cds != "") %>%
    select(allele, exon, cds)

write_tsv(out_df, paste0("./index_", locus, ".tsv"))
