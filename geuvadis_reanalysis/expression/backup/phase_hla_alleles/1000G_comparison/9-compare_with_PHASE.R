devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)


hap_diffs <- read_tsv("./hla_diffs_to_1000Ghaps.tsv")


kgp <- read_tsv("./hla_haps_phased.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(locus = sub("HLA-", "", locus)) 

kgp_w <- kgp %>%
    spread(locus, allele)

phase <- 
    "~/hlaexpression/geuvadis_reanalysis/phase_hla/phase_hla_haps_snps.tsv" %>%
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
   
write_tsv(concordant_set, "./concordant_set.tsv")
