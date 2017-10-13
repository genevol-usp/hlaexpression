devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- geuvadis_info %>% 
    filter(pop != "YRI", kgp_phase3 == 1) %>%
    pull(name)

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_id, gene_name)

pri_qtls <- 
    "../../../pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_best == 1L) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rank, variant = var_id)

imgt_qtls <- 
    "../../3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_best == 1L) %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(gene = gene_name, rank, variant = var_id)

qtls_df <- bind_rows(list(imgt = imgt_qtls, pri = pri_qtls), .id = "source")

writeLines(unique(c(qtls_df$variant)), "./imgt_and_pri.txt")

system("./subset_vcf.sh")

vcf <- read_tsv("./hla_eqtls.vcf", comment = "##") %>%
    select(ID, samples) %>%
    gather(subject, genotype, -1)

haps_df <- left_join(vcf, qtls_df, by = c("ID" = "variant")) %>%
    select(gene, source, rank, variant = ID, subject, genotype)

genotypes_df <- haps_df %>%
    separate(genotype, c("1", "2"), sep = "\\|") %>%
    gather(hap, allele, `1`:`2`) %>%
    arrange(subject, gene, source, rank, hap)

haps_f <- 
    genotypes_df %>%
    select(-variant) %>%
    spread(source, allele) %>%
    count(gene, rank, pri, imgt) %>%
    group_by(gene, rank) %>%
    mutate(f_hap = n/sum(n)) %>%
    ungroup() %>%
    select(gene, rank, allele_pri = pri, allele_imgt = imgt, f_hap)

alleles_f <- genotypes_df %>%
    count(gene, source, rank, allele) %>%
    group_by(gene, rank, source) %>%
    mutate(f = n/sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    spread(source, f) %>%
    rename(f_pri = pri, f_imgt = imgt)

ld_df <-
    left_join(haps_f, select(alleles_f, gene, rank, allele, f_pri), 
	      by = c("gene", "rank", "allele_pri" = "allele")) %>%
    left_join(select(alleles_f, gene, rank, allele, f_imgt), 
	      by = c("gene", "rank", "allele_imgt" = "allele")) %>%
    mutate(D = f_hap - f_pri*f_imgt,
	   rsq = D^2 / ((f_pri * (1 - f_pri)) * (f_imgt * (1 - f_imgt))))

write_tsv(ld_df, "./LD_imgt_pri_bestSNPs.tsv")

