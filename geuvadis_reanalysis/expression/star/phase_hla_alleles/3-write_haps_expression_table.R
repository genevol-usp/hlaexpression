devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

pops_df <- select(geuvadis_info, subject = name, pop)

allele_codes <- read_tsv("./PHASE/codes-phase.inp")

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    mutate(locus = sub("^HLA-", "", gene_name))

PHASE_calls <-
    read_phase("./PHASE/phase.out", loci = gencode_hla$locus) %>%
    rename(code = allele) %>%
    left_join(allele_codes, by = c("locus", "code")) %>%
    left_join(pops_df) %>%
    filter(pop != "YRI") %>%
    select(subject, locus, hap, hla_allele = allele)

kgp_calls <- 
    read_tsv("./1000G_comparison/hla_haps_mapped_to_1000G_phased.tsv") %>% 
    select(-S) %>%
    gather(locus, hla_allele, A:DRB1) %>% 
    mutate(hap = sub("^a", "", hap)) %>% 
    left_join(pops_df) %>%
    filter(pop != "YRI") %>%
    select(subject, locus, hap, hla_allele)

PHASE_haps <-
    PHASE_calls %>%
    spread(locus, hla_allele) %>%
    select(-hap) %>%
    distinct()
  
kgp_haps <- spread(kgp_calls, locus, hla_allele)

concordant_haps <- inner_join(kgp_haps, PHASE_haps)

concordant_class1 <- 
    inner_join(kgp_haps %>% select(subject:C), 
	       PHASE_haps %>% select(subject:C) %>% distinct())

concordant_class2 <- 
    inner_join(kgp_haps %>% select(subject, hap, DQA1:DRB1), 
	       PHASE_haps %>% select(subject, DQA1:DRB1) %>% distinct)

write_tsv(concordant_class1, "./data/concordant_haps_classI.tsv")
write_tsv(concordant_class2, "./data/concordant_haps_classII.tsv")
write_tsv(concordant_haps, "./data/concordant_haps_classIandII.tsv")

expression_df <-
    read_tsv("../quantifications_2/processed_quant.tsv") %>%
    mutate(subject = convert_ena_ids(subject),
	   locus = sub("^HLA-", "", locus),
	   allele = sub("IMGT_", "", allele)) %>%
    filter(locus %in% gencode_hla$locus, subject %in% kgp_haps$subject) %>%
    select(subject, locus, allele, tpm) %>%
    arrange(subject, locus)

eqtl_df <-
    read_tsv("./best_eqtls.tsv") %>%
    left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(locus, rank, variant = var_id)

eqtl_info <- 
    read_tsv("./best_eqtl_snps.vcf", comment = "##") %>%
    select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -ID) %>%
    inner_join(eqtl_df, by = c("ID" = "variant")) %>%
    select(subject, locus, rank, variant = ID, genotype) %>%
    separate(genotype, c("1", "2"), sep = "\\|") %>%
    gather(hap, allele, `1`:`2`) %>%
    arrange(subject, locus, rank, variant, hap)

kgp_calls_and_expr <- 
    kgp_calls %>%
    separate_rows(hla_allele, sep = "/IMGT_") %>%
    mutate(hla_allele = sub("IMGT_", "", hla_allele)) %>%
    group_by(subject, locus) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    arrange(subject, locus, hap) %>%
    left_join(expression_df, 
	      by = c("subject", "locus", "hla_allele" = "allele")) %>%
    distinct() %>%
    left_join(eqtl_info, by = c("subject", "locus", "hap")) %>%
    arrange(subject, locus, hap, hla_allele, rank) %>%
    rename(variant_allele = allele)

PHASE_calls_and_expr <-
    PHASE_calls %>%
    mutate(hla_allele = sub("IMGT_", "", hla_allele)) %>%
    left_join(expression_df, 
	      by = c("subject", "locus", "hla_allele" = "allele")) %>%
    distinct() 

write_tsv(kgp_calls_and_expr, "./data/1000G_haps_expression_snps.tsv")
write_tsv(PHASE_calls_and_expr, "./data/PHASE_haps_expression_snps.tsv")
