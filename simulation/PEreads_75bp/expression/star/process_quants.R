devtools::load_all("~/hlaseqlib")
library(tidyverse)

make_genot_calls_df <- function(typings_df) {

    typings_df %>%
	mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3))
}

quant_round <- commandArgs(TRUE)[1]

samples <- sprintf("sample_%02d", 1:50)
    
quant_files <- 
    paste0("./quantifications_", quant_round) %>%
    file.path(samples, "quant_imgt.sf") %>%
    setNames(samples)

missing_files <- quant_files[!file.exists(quant_files)]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

hla_genes <- sort(gencode_hla$gene_name)

if (quant_round == 1L | quant_round == 2L) {

    quants <- quant_files %>%
	plyr::ldply(read_star_imgt_quants, .id = "subject") %>%
	mutate(subject = as.character(subject))
    
    goldstd_genos <- read_tsv("../../data/genos.tsv")

    if (quant_round == 1L) {
    
	thresholds <- as.list(seq(0, .25, .05))
	names(thresholds) <- seq(0, .25, .05)

	typings <- 
	    plyr::ldply(thresholds, 
			function(th) hla_genotype_dt(quants, th),
			.id = "th")

	calls <- 
	    typings %>%
	    filter(locus %in% hla_genes) %>%
	    select(th, subject, locus, allele) %>%
	    make_genot_calls_df()

	accuracies <-
	    calls %>%
	    split(.$th) %>%
	    plyr::ldply(function(df) calc_genotyping_accuracy(df, goldstd_genos),
			.id = "th") %>%
	    group_by(th) %>%
	    mutate(th_average = mean(accuracy)) %>%
	    ungroup()
	
	best_th <- accuracies %>%
	    slice(which.max(th_average)) %>%
	    pull(th)

	write_tsv(accuracies, "./genotyping_accuracies_1.tsv")

	out_df <- filter(typings, th == best_th) %>% select(-th)

  } else if (quant_round == 2L) {

    out_df <- hla_genotype_dt(quants, th = 0) %>%
	hla_apply_zigosity_threshold(th = 0.25)

    calls <- out_df %>%
	filter(locus %in% hla_genes) %>%
	select(subject, locus, allele) %>%
	make_genot_calls_df()

    accuracies <- calc_genotyping_accuracy(calls, goldstd_genos)

    write_tsv(accuracies, "./genotyping_accuracies_2.tsv")  }

} else if (quant_round == "PRI")  {

  out_df <- quant_files %>%
      plyr::ldply(read_star_pri_quants, .id = "subject")
}

out_df %>%
    write_tsv(paste0("./quantifications_", quant_round, "/processed_quant.tsv"))

