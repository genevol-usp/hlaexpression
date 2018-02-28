library(Biostrings)
library(tidyverse)
devtools::load_all("/home/vitor/hlaseqlib")
  
set.seed(2)
sample_df <- geuvadis_info %>%
    filter(kgp_phase3 == 1L, pop != "YRI") %>%
    select(sample_id = name, subject = ena_id) %>%
    filter(sample_id %in% sample(sample_id, 50)) %>%
    arrange(sample_id) %>%
    mutate(code = sprintf("sample_%02d", 1:50)) %>%
    select(subject, code)

abundances <- 
    file.path("~/hlaexpression/geuvadis_reanalysis/expression/star/supplemented/quantifications_2",
	      sample_df$subject, "quant.sf") %>%
    setNames(sample_df$subject) %>%
    map_df(read_tsv, .id = "subject") %>%
    left_join(sample_df, by = "subject") %>%
    select(subject = code, Name, Length, EffectiveLength, TPM, NumReads)

genos <- abundances %>%
    filter(grepl("IMGT", Name)) %>%
    mutate(locus = imgt_to_gname(Name)) %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele = Name) %>%
    mutate(allele = sub("IMGT_", "", allele),
	   allele = hla_trimnames(allele, 3)) %>%
    group_by(subject, locus) %>%
    mutate(n = ifelse(n() == 1L, "1_1", "1")) %>%
    ungroup() %>%
    separate_rows(n, sep = "_") %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus, allele)

write_tsv(genos, "./genos.tsv")

index <- readDNAStringSet("../../../imgt_index_v2/gencode.v25.PRI.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

tx <- intersect(names(index), unique(abundances$Name))

abundances_df_intersect <- abundances %>%
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
    mutate_at(vars(sample_01:sample_50), ~ifelse(is.na(.), 0L, .))

index <- index[phenotypes$Name]
writeXStringSet(index, "./index_simulation.fa")

phenotypes %>%
    select(-Name) %>%
    write_tsv("./phenotypes.tsv")
