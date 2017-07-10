library(data.table)
devtools::load_all("~/hlatools")
devtools::load_all("/home/vitor/genomicRutils/")

process_quant_imgt <- function(quant_files) {

  readDT <- function(x)
    fread(x, select = c(1, 4, 5))[grepl("^IMGT", target_id)]

  ldt <- parallel::mclapply(quant_files, readDT, mc.cores = 48)

  rbindlist(ldt, idcol = "subject")[, .(subject, target_id, est_counts, tpm)]
}

samples <- readLines("./nci_samples.txt")

quant_round <- commandArgs(TRUE)[1]
quant_dir <- paste0("./quantifications_", quant_round)
quant_files <- file.path(quant_dir, samples, "abundance.tsv")
names(quant_files) <- samples
  
if (quant_round == 1 | quant_round == 2) {
  
  imgt_quants <- process_quant_imgt(quant_files)
  
  hla_names <- imgt_quants[, .(target_id = unique(target_id))]
  hla_names[, `:=`(locus = imgt_to_gname(target_id),
		   gene_id = imgt_to_gid(target_id))]

  imgt_quants <- 
    imgt_quants[hla_names, on = .(target_id)
	      ][, .(subject, gene_id, locus, allele = target_id, est_counts, tpm)]
  
  if (quant_round == 1) {

    thresholds <- as.list(seq(0, .25, .05))
    names(thresholds) <- seq(0, .25, .05)

    typings <-
      parallel::mclapply(thresholds, 
			 function(th) hla_genotype_dt(imgt_quants, th),
			 mc.cores = length(thresholds))
      
    typings_dt <- rbindlist(typings, idcol = "th")

    genos_dt <- typings_dt[locus %in% paste0("HLA-", c("A", "B", "C", "DQB1", "DRB1")), 
			   .(th, subject, locus, allele)]
    genos_dt[, allele := hla_trimnames(gsub("IMGT_|_s\\d", "", allele), 3)]

    nci_genos <- fread("../data/nci_expression.tsv", select = c(1, 2, 4))
   
    accuracies_list <- 
      parallel::mclapply(split(genos_dt, genos_dt$th),
			 function(i) calc_genotyping_accuracy(i, nci_genos),
			 mc.cores = length(thresholds))

    accuracies_dt <- rbindlist(accuracies_list, idcol = "th")
    accuracies_dt[, th_average := mean(accuracy), by = th] 

    best_th <- accuracies_dt[which.max(th_average), th]

    accuracies_dt[, `:=`(accuracy = round(accuracy, 2), 
			 th_average = round(th_average, 2))]

    fwrite(accuracies_dt, "./genotyping_accuracies.tsv", 
	   sep = "\t", quote = FALSE)

    out_dt <- typings_dt[th == best_th]
    out_dt[, th := NULL]

    best_genos <- genos_dt[th == best_th]
    best_genos[, th := NULL]

    best_calls <- calc_genotyping_accuracy(best_genos, nci_genos, FALSE)
    setnames(best_calls, c("allele.x", "allele.y", "correct"), 
	     c("allele_rnaseq", "allele_nci", "concordant"))
    fwrite(best_calls, "./allele_calls.tsv", sep = "\t", quote = FALSE)

  } else if (quant_round == 2) {

  out_dt <- hla_genotype_dt(imgt_quants, th = 0)

  }

} else {

  stop("wrong dir")
 
}

fwrite(out_dt, file.path(quant_dir, "processed_quant.tsv"), 
       sep = "\t", quote = FALSE)
