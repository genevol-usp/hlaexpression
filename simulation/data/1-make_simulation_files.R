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

processed_imgt <- 
    "~/hlaexpression/geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    inner_join(sample_df, by = "subject") %>%
    select(subject = code, Name = allele, NumReads = est_counts)

abundances <- 
    file.path("~/hlaexpression/geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications",
	      sample_df$subject, "quant.sf") %>%
    setNames(sample_df$subject) %>%
    map_df(read_tsv, .id = "subject") %>%
    left_join(sample_df, by = "subject") %>%
    select(subject = code, Name, NumReads) %>%
    filter(!grepl("IMGT_", Name)) %>%
    bind_rows(distinct(processed_imgt)) %>%
    arrange(subject, Name)

genos <- processed_imgt %>%
    mutate(locus = imgt_to_gname(Name),
	   allele = gsub("IMGT_", "", Name)) %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele) %>%
    arrange(subject, locus)
    
write_tsv(genos, "./genos.tsv")

index <- readDNAStringSet("../../imgt_index/gencode.v25.PRI.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

tx <- intersect(names(index), unique(abundances$Name))

abundances_df_intersect <- abundances %>%
    filter(Name %in% tx) %>%
    group_by(subject) %>%
    mutate(TrueCounts = as.integer(round(NumReads/sum(NumReads) * 3e7))) %>%
    ungroup()

phenotype_counts <- abundances_df_intersect %>%
    select(subject, Name, TrueCounts)

write_tsv(phenotype_counts, "./phenotypes_trueCounts.tsv")

phenotypes <- phenotype_counts %>%
    spread(subject, TrueCounts) %>%
    mutate_at(vars(sample_01:sample_50), ~ifelse(is.na(.), 0L, .))

index <- index[phenotypes$Name]
writeXStringSet(index, "./index_simulation.fa")

phenotypes %>%
    select(-Name) %>%
    write_tsv("./phenotypes.tsv")
