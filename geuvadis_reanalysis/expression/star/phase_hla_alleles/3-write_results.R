devtools::load_all("~/hlaseqlib")
library(tidyverse)

expression_df <- 
    read_tsv("../supplemented/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    left_join(geuvadis_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    select(subject = name, locus, allele, tpm) %>%
    arrange(subject, locus)

kgp_calls <- read_tsv("./1000G_comparison/hla_haps_phased.tsv") 
  

# I removed this DQB1 allele bc:
# This is actually a pair of alleles which share the same 3 fields but don't
# have the same CDS, because in insertion in on of them.
# Because in the haplotype assignment I only consider positions in the ref
# genome, they are equal. So the individual having these alleles is homozygote
# in this respect.
# These alleles are probably associated to typing errors. True heterozygotes for
# one of the alleles are called as heterozygotes with the 2 alleles, with the
# other true allele missing. This happens because these alleles are so similar
# that the reads of one map to the other.

haps_exp_df <- kgp_calls %>%
    filter(!grepl("/", allele), !grepl("DQB1\\*05:03:01", allele)) %>%
    left_join(expression_df, by = c("subject", "locus", "allele")) %>% 
    distinct(subject, locus, hap, allele, .keep_all = TRUE) %>%
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
