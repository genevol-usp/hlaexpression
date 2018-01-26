devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_ref <- read_tsv("~/hlaexpression/imgt_index/hla_ref_alleles.tsv") %>%
    mutate(subject = "ERR188021",
	   allele = hla_trimnames(allele, 3)) %>%
    select(subject, locus, allele)

gold_std <- bind_rows(hla_ref, hla_ref) %>% arrange(locus)

imgt_quants <- read_tsv("./quantifications_1/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name, 
	   est_counts = NumReads, tpm = TPM)

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- 
  plyr::ldply(thresholds, 
	      function(th) hla_genotype_dt(imgt_quants, th),
	      .id = "th")

calls <- typings %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(th, subject, locus, allele) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3))

accuracies <- calls %>%
    split(.$th) %>%
    plyr::ldply(function(df) calc_genotyping_accuracy(df, gold_std),
		.id = "th") %>%
    group_by(th) %>%
    mutate(th_average = mean(accuracy)) %>%
    ungroup()
  
write_tsv(accuracies, "./genotyping_accuracies_1.tsv")

best_th <- accuracies %>%
  slice(which.max(th_average)) %>%
  pull(th) %>%
  as.character()

out_df <- filter(typings, th == best_th) %>% select(-th)

write_tsv(out_df, "./quantifications_1/processed_imgt_quants.tsv")
