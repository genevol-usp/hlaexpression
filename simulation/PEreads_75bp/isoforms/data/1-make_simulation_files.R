library(tidyverse)
library(Biostrings)
devtools::load_all("/home/vitor/hlaseqlib")

annot <-
    read_tsv("~/gencode_data/gencode.v25.annotation.gtf", 
             comment = "##", col_names = FALSE,
             col_types = "c-cii-c-c", progress = FALSE) %>%
    filter(X3 == "transcript") %>% 
    select(X9) %>%
    mutate(i = seq_len(nrow(.)), X9 = stringr::str_split(X9, "; ")) %>%
    unnest() %>%
    filter(grepl("^gene_name|^gene_id|^transcript_id|^transcript_type", X9)) %>%
    separate(X9, c("tag", "id"), " ") %>%
    mutate(id = gsub("\"", "", id)) %>%
    spread(tag, id) %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, transcript_id, transcript_type)

abundances <- 
    "../../../../geuvadis_reanalysis/expression/star/pri/quantifications/ERR188021/quant.sf" %>%
    read_tsv()

transcripts_to_rm <- abundances %>%
    inner_join(annot, by = c("Name" = "transcript_id")) %>%
    filter(TPM < 10 | transcript_type != "protein_coding") %>%
    pull(Name)

abundances <- abundances %>%
    filter(! Name %in% transcripts_to_rm)

index <- readDNAStringSet("~/gencode_data/gencode.v25.PRI.transcripts.fa")
index <- index[width(index) >= 75]

tx <- intersect(names(index), abundances$Name)

index <- index[tx]
writeXStringSet(index, "./polyester_index.fa")

abundances <- abundances %>%
    filter(Name %in% tx) %>%
    mutate(TrueCounts = as.integer(round(NumReads/sum(NumReads) * 3e7)),
	   TrueTPM = counts_to_tpm(TrueCounts, EffectiveLength))

phenotypes_tpm_counts <- abundances %>%
    select(Name, TrueCounts, TrueTPM)

write_tsv(phenotypes_tpm_counts, "./phenotypes_counts_tpm.tsv")

phenotypes <- phenotypes_tpm_counts %>%
    select(sample_01 = TrueCounts)

write_tsv(phenotypes, "./phenotypes.tsv")
