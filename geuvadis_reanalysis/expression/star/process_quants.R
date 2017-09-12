devtools::load_all("~/hlaseqlib")
library(tidyverse)

doMC::registerDoMC(25)

make_genot_calls_df <- function(typings_df) {
    
    typings_df %>%
        mutate(subject = convert_ena_ids(as.character(subject)),
	       locus = sub("^HLA-", "", locus),
	       allele = hla_trimnames(gsub("IMGT_", "", allele), 3))
}

read_imgt_quants <- function(f) {

    read_tsv(f, col_types = "c--dd", progress = FALSE) %>% 
	filter(grepl("^IMGT_", Name)) %>%
	mutate(locus = imgt_to_gname(Name), gene_id = gname_to_gid(locus)) %>%
	select(locus, gene_id, allele = Name, est_counts = NumReads, tpm = TPM)
}

read_pri_quants <- function(f) {
    read_tsv(f, col_types = "c--dd", progress = FALSE) %>%
    separate(Name, c("tx_id", "gene_id", "dummy1", "dummy2", "tx_name", 
		     "gene_name", "dummy3", "dummy4", "dummy5"), 
	     sep = "\\|") %>%
    select(tx_id, gene_name, est_counts = NumReads, tpm = TPM) %>%
    filter(gene_name %in% hla_genes) %>%
    group_by(gene_name) %>%
    summarize_at(vars(tpm, est_counts), sum) %>%
    ungroup() %>%
    rename(locus = gene_name)
}

quant_round <- as.integer(commandArgs(TRUE)[1])

samples <- 
    geuvadis_info %>% 
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

quant_files <- paste0("./quantifications_", quant_round) %>%
      file.path(samples, "quant.sf") %>%
      setNames(samples)

missing_files <- quant_files[!file.exists(quant_files)]

if (length(missing_files) == 0L) {
    stop(paste("missing files:", missing_files))
}

hla_genes <- paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))

if (quant_round == 1L || quant_round == 2L) {

    quants <- 
	quant_files %>%
	plyr::ldply(read_imgt_quants, .id = "subject", .parallel = TRUE)
    
    goldstd_genos <- mutate(pag, allele = hla_trimnames(allele, 3))
    
    if (quant_round == 1L) {
    
        thresholds <- as.list(seq(0, .25, .05))
        names(thresholds) <- seq(0, .25, .05)

        typings <- 
	  plyr::ldply(thresholds, function(th) hla_genotype_dt(quants, th),
		      .id = "th", .parallel = TRUE)

	calls <- 
	    typings %>%
	    filter(locus %in% hla_genes) %>%
	    select(th, subject, locus, allele) %>%
	    make_genot_calls_df() 

        accuracies <-
            calls %>%
            split(.$th) %>%
            plyr::ldply(function(df) calc_genotyping_accuracy(df, genos),
			.id = "th", .parallel = TRUE) %>%
            group_by(th) %>%
            mutate(th_average = mean(accuracy)) %>%
            ungroup()
	  
      write_tsv(accuracies, "./genotyping_accuracies_1.tsv")

      best_th <- accuracies %>%
          slice(which.max(th_average)) %>%
          pull(th) %>%
          as.character()

      out_df <- filter(typings, th == best_th) %>% select(-th)

    } else if (quant_round == 2L) {

        out_df <- hla_genotype_dt(quants, th = 0) %>%
	    hla_apply_zigosity_threshold(th = 0.25)
    
	calls <- 
	    out_df %>%
	    filter(locus %in% hla_genes) %>%
	    select(th, subject, locus, allele) %>%
	    make_genot_calls_df()
	
        accuracies <- calc_genotyping_accuracy(calls, goldstd_genos)
      
	write_tsv(accuracies, "./genotyping_accuracies_2.tsv")
    }

} else if (quant_round == "PRI") {

    out_df <- 
	quant_files %>%
	plyr::ldply(read_pri_quants, .id = "subject", .parallel = TRUE)
}

out_df %>%
    write_tsv(paste0("./quantifications_", quant_round, "/processed_quant.tsv"))

