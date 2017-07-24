devtools::load_all("~/hlaseqlib")
library(tidyverse)

quant_round <- commandArgs(TRUE)[1]

genos <- read_tsv("../../data/genos.tsv")

samples <- sprintf("sample_%02d", 1:50)

doMC::registerDoMC(25)

quants <- 
  file.path(paste0("./quantifications_", quant_round), samples, "quant.sf") %>%
  setNames(samples) %>%
  plyr::ldply(. %>% read_tsv(col_types = "c--dd", progress = FALSE) %>% 
	      filter(grepl("^IMGT_", Name)) %>%
	      mutate(locus = imgt_to_gname(Name),
		     gene_id = gname_to_gid(locus)) %>%
	      select(locus, gene_id, allele = Name, est_counts = NumReads, tpm = TPM),
	    .id = "subject", .parallel = TRUE)

if (quant_round == 1L) {
  
  thresholds <- as.list(seq(0, .25, .05))
  names(thresholds) <- seq(0, .25, .05)

  typings <- plyr::ldply(thresholds, function(th) hla_genotype_dt(quants, th),
			 .id = "th", .parallel = TRUE)

  calls <- typings %>%
    filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
    select(th, subject, locus, allele) %>%
    mutate(subject = as.character(subject),
	   allele = hla_trimnames(gsub("IMGT_", "", allele), 3))

  accuracies <-
    calls %>%
    split(.$th) %>%
    plyr::ldply(function(df) calc_genotyping_accuracy(df, genos),
		.id = "th", .parallel = TRUE) %>%
    group_by(th) %>%
    mutate(th_average = mean(accuracy)) %>%
    ungroup()

  best_th <- accuracies %>%
    slice(which.max(th_average)) %>%
    pull(th)

  accuracies <- accuracies %>%
    mutate(accuracy = round(accuracy, 2),
	   th_average = round(th_average, 3))

  write_tsv(accuracies, "./genotyping_accuracies.tsv")

  out_df <- typings %>%
    filter(th == best_th) %>%
    select(-th)

} else if (quant_round == 2L) {

  out_df <- hla_genotype_dt(quants, th = 0)
}

write_tsv(out_df, paste0("./quantifications_", quant_round, "/processed_quant.tsv"))
