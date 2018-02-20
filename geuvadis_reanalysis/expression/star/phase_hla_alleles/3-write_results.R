devtools::load_all("~/hlaseqlib")
library(tidyverse)

expression_df <- 
    read_tsv("../supplemented/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    left_join(geuvadis_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = gsub("IMGT_", "", allele),
	   allele_3F = hla_trimnames(allele, 3)) %>%
    select(subject = name, locus, allele, allele_3F, tpm) %>%
    arrange(subject, locus)

kgp_calls <- read_tsv("./1000G_comparison/hla_haps_phased.tsv") 
  

# DQB1*05:03:01:01 and DQB1*05:03:01:03:
# Pair of alleles which share the same 3 fields but don't have the same CDS due
# to an insertion in on of them. Because in the haplotype assignment I only
# consider positions in the ref genome, they are equal. So the individual having
# these alleles is homozygote in this respect. These alleles are associated
# with typing errors. True heterozygotes for one of these alleles are called as
# heterozygotes with the 2 alleles, with the other true allele missing. This
# happens because these alleles are so similar that the reads of one map to the
# other.

haps_exp_df <- kgp_calls %>%
    separate_rows(allele, sep = "/") %>%
    left_join(expression_df, by = c("subject", "locus", "allele")) %>%
    left_join(expression_df, by = c("subject", "locus", "allele" = "allele_3F")) %>%
    mutate(tpm = ifelse(is.na(tpm.x), tpm.y, tpm.x)) %>%
    distinct(subject, locus, hap, allele.y, tpm) %>%
    group_by(subject, locus, hap) %>%
    summarise(allele = paste(allele.y, collapse = "/"),
	      tpm = paste(tpm, collapse = "/")) %>%
    ungroup() %>%
    rename(hla_allele = allele)

eqtl_df <- read_tsv("./best_eqtl.tsv") %>%
    left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(locus = gene_name, rsid = var_id)

eqtl_info <- read_tsv("./best_eqtl_snps.vcf", comment = "##") %>%
    select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -ID) %>%
    left_join(eqtl_df, by = c("ID" = "rsid")) %>%
    select(subject, locus, rsid = ID, genotype) %>%
    separate(genotype, c("1", "2"), sep = "\\|") %>%
    gather(hap, allele, `1`, `2`) %>%
    mutate(hap = as.integer(hap)) %>% 
    arrange(subject, locus, rsid, hap)

haps_exp_snp_df <- 
    left_join(haps_exp_df, eqtl_info, by = c("subject", "locus", "hap")) %>%
    select(subject, locus, hap, rsid, allele, hla_allele, tpm)

write_tsv(haps_exp_snp_df, "./data/1000G_haps_expression_rsid.tsv")
