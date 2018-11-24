devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

hap_diffs <- read_tsv("./kgp_comparison/hla_diffs_to_1000Ghaps.tsv")

kgp <- read_tsv("./kgp_comparison/hla_haps_phased.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(locus = sub("HLA-", "", locus)) 

kgp_w <- kgp %>%
    spread(locus, allele)

phase <- 
    "./phase/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>%
    distinct(subject, locus, hap, allele = allele_gene, unsure = uncertain_gene) %>%
    mutate(locus = sub("HLA-", "", locus)) 

phase_w <- phase %>%
    select(-unsure) %>%
    spread(locus, allele) %>%
    select(-hap)

concordant_set <- inner_join(kgp_w, phase_w) %>%
    select(subject, hap) %>%
    inner_join(kgp) %>%
    mutate(locus = factor(paste0("HLA-", locus), levels = gencode_hla$gene_name)) %>%
    arrange(subject, locus, hap)
  
phase_snps <- 
    "./phase/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>%
    filter(rank == 0) %>%
    select(subject, locus, allele = allele_gene, qtl_rsid = rsid, 
	   qtl_allele = allele_snp)

hlapers <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele)) %>%
    select(subject, locus, allele, tpm)

out <- 
    left_join(concordant_set, phase_snps, by = c("subject", "locus", "allele")) %>%
    distinct(subject, hap, locus, .keep_all = TRUE) %>%
    left_join(hlapers, by = c("subject", "locus", "allele")) %>%
    distinct(subject, hap, locus, .keep_all = TRUE)

write_tsv(out, "./concordant_set.tsv")
