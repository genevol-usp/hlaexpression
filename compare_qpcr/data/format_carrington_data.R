library(tidyverse)
library(readxl)

pcr_data <- read_excel("./HLA-A, -B, -C expression levels to Vitor.xlsx") %>%
    select(subject = 1, A = 2, B = 3, C = 4, 
           allele_A1 = 6, allele_A2 = 7, allele_B1 = 8, allele_B2 = 9,
           allele_C1 = 10, allele_C2 = 11)
        
pcr_rna <- select(pcr_data, subject, A:C) %>%
    gather(locus, rna, -1) %>%
    mutate(rna = as.numeric(rna)) %>%
    filter(!is.na(rna)) %>%
    arrange(subject, locus)

pcr_genos <- select(pcr_data, subject, allele_A1:allele_C2) %>%
    gather(locus, allele, -1) %>%
    mutate(locus = sub("^allele_([ABC])[1-2]$", "\\1", locus),
           lineage = sub("^([^:]+).*$", "\\1", allele),
           lineage = paste(locus, lineage, sep = "*")) %>%
    select(subject, locus, lineage)

pcr <- inner_join(pcr_genos, pcr_rna, by = c("subject", "locus")) %>%
    select(subject, locus, lineage, rna) %>%
    arrange(subject, locus, lineage)

write_tsv(pcr, "./hla_carrington_data.tsv")

pcr_filtered <- pcr %>% 
    group_by(lineage) %>%
    filter(n_distinct(subject) >= 10) %>%
    ungroup() %>% 
    arrange(subject, locus, lineage) %>%
    select(subject, locus, lineage, rna)

write_tsv(pcr_filtered, "./hla_carrington_data_filteredN10.tsv")
