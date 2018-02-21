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

#kgp_haps <- kgp_calls %>%
#    mutate(locus = sub("HLA-", "", locus)) %>%
#    spread(locus, allele)
#
#kgp_diff <- read_tsv("./1000G_comparison/hla_diffs_to_1000Ghaps.tsv")
#
#allele_codes <- read_tsv("./PHASE/codes-phase.inp") %>%
#    mutate(allele = sub("IMGT_", "", allele))
#
#phase_calls <-
#    read_phase("./PHASE/phase.out", 
#	       loci = sub("HLA-", "", gencode_hla$gene_name)) %>%
#    rename(code = allele) %>%
#    left_join(allele_codes, by = c("locus", "code")) %>%
#    select(subject, locus, hap, hla_allele = allele, uncertain)
#
#phase_haps <- phase_calls %>%
#    select(-uncertain) %>%
#    spread(locus, hla_allele) %>%
#    select(-hap)
#
#kgp_class1 <- select(kgp_haps, subject:C)
#phase_class1 <- select(phase_haps, subject:C)
#
#diff_subjects <- anti_join(kgp_class1, phase_class1) %>%
#    filter(!grepl("/", A), !grepl("/", B), !grepl("/", C)) %>%
#    distinct(subject)
#
#s <- diff_subjects$subject[2] 
#kgp_class1 %>% filter(subject == s)
#phase_class1 %>% filter(subject == s)
#kgp_diff %>% filter(subject == s, locus %in% paste0("HLA-", c("A", "B", "C")))



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
    mutate(allele = ifelse(is.na(allele.y), allele, allele.y),
	   tpm = ifelse(is.na(tpm.x), tpm.y, tpm.x)) %>%
    distinct(subject, locus, hap, allele, tpm) %>%
    group_by(subject, locus, hap) %>%
    summarise(allele = paste(allele, collapse = "/"),
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
