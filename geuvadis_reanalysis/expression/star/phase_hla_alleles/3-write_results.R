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
   
haps_exp_df <- kgp_calls %>%
    filter(!grepl("/", allele)) %>%
    left_join(expression_df, by = c("subject", "locus", "allele")) %>%
    distinct(subject, locus, hap, allele, .keep_all = TRUE) %>%
    left_join(expression_df, by = c("subject", "locus", "allele" = "allele_3F")) %>%
    distinct(subject, locus, hap, allele.y, .keep_all = TRUE) %>%
    mutate(tpm = ifelse(is.na(tpm.x), tpm.y, tpm.x)) %>%
    select(subject, locus, hap, hla_allele = allele.y, tpm)

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

write_tsv(out_df, "./data/1000G_haps_expression_rsid.tsv")
