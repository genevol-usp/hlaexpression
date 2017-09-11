devtools::load_all("~/hlaseqlib")
library(tidyverse)

quant_round <- as.integer(commandArgs(TRUE)[1])

samples <- 
  geuvadis_info %>% 
  filter(kgp_phase3 == 1L & pop != "YRI") %>%
  pull(ena_id)

doMC::registerDoMC(25)

if (quant_round == 1L || quant_round == 2L) {

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
    
    genos <- mutate(pag, allele = hla_trimnames(allele, 3))
    
    thresholds <- as.list(seq(0, .25, .05))
    names(thresholds) <- seq(0, .25, .05)

    typings <- plyr::ldply(thresholds, function(th) hla_genotype_dt(quants, th),
			   .id = "th", .parallel = TRUE)

    calls <- typings %>%
      filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
      select(th, subject, locus, allele) %>%
      mutate(subject = convert_ena_ids(as.character(subject)),
	     locus = sub("^HLA-", "", locus),
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
      pull(th) %>%
      as.character()

    accuracies <- accuracies %>%
      mutate(accuracy = round(accuracy, 2),
	     th_average = round(th_average, 3))

    write_tsv(accuracies, "./genotyping_accuracies.tsv")

    out_df <- typings %>%
      filter(th == best_th) %>%
      select(-th)

  } else if (quant_round == 2L) {

    out_df <- hla_genotype_dt(quants, th = 0) %>%
      mutate(subject = as.character(subject)) %>%
      group_by(subject, locus) %>%
      mutate(i = as.integer(any(tpm/max(tpm) <= 0.25))) %>%
      ungroup()

    out_df_0 <- out_df %>% filter(i == 0L) %>% select(-i)

    out_df_1 <- 
      out_df %>% 
      filter(i == 1L) %>%
      group_by(subject, gene_id) %>%
      mutate(m = as.integer(tpm == max(tpm)), 
             est_counts = sum(est_counts)/2, tpm = sum(tpm)/2) %>%
      ungroup() %>%
      filter(m == 1L) %>%
      select(-i, -m)

    out_df <- bind_rows(out_df_0, out_df_1, out_df_1) %>%
      mutate(subject = as.character(subject)) %>%
      arrange(subject, locus, allele)
    
   # genos <- mutate(pag, allele = hla_trimnames(allele, 3))
   # 
   # calls <- out_df %>%
   #   filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
   #   select(subject, locus, allele) %>%
   #   mutate(subject = convert_ena_ids(as.character(subject)),
   #          locus = sub("^HLA-", "", locus),
   #          allele = hla_trimnames(gsub("IMGT_", "", allele), 3))

   # accuracy <- calc_genotyping_accuracy(calls, genos) %>%
   #   mutate(accuracy = round(accuracy, 3))

} else if (quant_round == "CHR") {

  out_df <- 
    file.path(paste0("./quantifications_", quant_round), samples, "quant.sf") %>%
    setNames(samples) %>%
    plyr::ldply(. %>%
		read_tsv(col_types = "c--dd", progress = FALSE) %>%
		separate(Name, c("tx_id", "gene_id", "dummy1", "dummy2", "tx_name", "gene_name", 
				 "dummy3", "dummy4", "dummy5"), sep = "\\|") %>%
		select(tx_id, gene_name, est_counts = NumReads, tpm = TPM) %>%
		filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
		group_by(gene_name) %>%
		summarize_at(vars(tpm, est_counts), sum) %>%
		ungroup() %>%
		rename(locus = gene_name),
	      .id = "subject", .parallel = TRUE)
}

if (!all(samples %in% out_df$subject)) {
  
  miss <- samples[! samples %in% out_df$subject]
  stop(paste("missing samples:", miss))
}

write_tsv(out_df, paste0("./quantifications_", quant_round, "/processed_quant.tsv"))
