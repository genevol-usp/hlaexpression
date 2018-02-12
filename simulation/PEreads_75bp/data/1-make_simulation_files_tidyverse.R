library(Biostrings)
library(tidyverse)
devtools::load_all("/home/vitor/hlaseqlib")
  
set.seed(2)
sample_df <- geuvadis_info %>%
    filter(kgp_phase3 == 1L, pop != "YRI") %>%
    select(sample_id = name, subject = ena_id) %>%
    filter(sample_id %in% sample(sample_id, 50)) %>%
    arrange(sample_id) %>%
    mutate(code = sprintf("sample_%02d", 1:50))

abundances <- 
    file.path("~/hlaexpression/geuvadis_reanalysis/expression/star/supplemented/quantifications_2",
	      sample_df$subject, "quant.sf") %>%
    setNames(sample_df$subject) %>%
    map_df(read_tsv, .id = "subject")

imgt_lens <- abundances %>%
    filter(grepl("IMGT", Name)) %>%
    select(subject, allele = Name, EffectiveLength) 

abundances_no_imgt <- abundances %>%
    filter(!grepl("IMGT", Name)) %>%
    left_join(sample_df, "subject") %>%
    select(subject = code, Name, EffectiveLength, NumReads) 
    
abundances_imgt <- 
    "~/hlaexpression/geuvadis_reanalysis/expression/star/supplemented/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    inner_join(sample_df, by = "subject") %>%
    left_join(imgt_lens, by = c("subject", "allele")) %>%
    select(subject = code, locus, Name = allele, EffectiveLength, NumReads = est_counts)

genos <- abundances_imgt %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele = Name) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", sub("^([^/]+).*$", "\\1", allele)), 3)) %>%
    arrange(subject, locus, allele)

write_tsv(genos, "./genos.tsv")

abundances_df <-
    select(abundances_imgt, -locus) %>%
    group_by(subject, Name) %>%
    summarise(EffectiveLength = unique(EffectiveLength),
	      NumReads = sum(NumReads)) %>%
    ungroup() %>%
    bind_rows(abundances_no_imgt, .) %>%
    arrange(subject, Name)

index <- readDNAStringSet("../../../imgt_index_v2/gencode.v25.PRI.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

tx <- intersect(names(index), unique(abundances_df$Name))

index <- index[tx]
writeXStringSet(index, "./polyester_index.fa")

abundances_df_intersect <- abundances_df %>%
    filter(Name %in% tx) %>%
    group_by(subject) %>%
    mutate(TrueCounts = as.integer(round(NumReads/sum(NumReads) * 3e7)),
	   TrueTPM = counts_to_tpm(TrueCounts, EffectiveLength)) %>%
    ungroup()

phenotypes_tpm_counts <- abundances_df_intersect %>%
    select(subject, Name, TrueCounts, TrueTPM)

write_tsv(phenotypes_tpm_counts, "./phenotypes_counts_tpm.tsv")

phenotypes <- phenotypes_tpm_counts %>%
    select(-TrueTPM) %>%
    spread(subject, TrueCounts) %>%
    mutate_at(vars(sample_01:sample_50), ~ifelse(is.na(.), 0L, .)) %>%
    select(-Name)

write_tsv(phenotypes, "./phenotypes.tsv")

phenotypes_old <- read_tsv("./phenotypes.tsv")
