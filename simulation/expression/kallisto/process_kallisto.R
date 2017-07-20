library(data.table)
devtools::load_all("~/hlatools")
devtools::load_all("/home/vitor/genomicRutils/")

process_quant_imgt <- function(quant_files) {

  readDT <- function(x)
    fread(x, select = c(1, 4, 5))[grepl("^IMGT", target_id)]

  ldt <- parallel::mclapply(quant_files, readDT, mc.cores = 50)

  rbindlist(ldt, idcol = "subject")[, .(subject, target_id, est_counts, tpm)]
}

process_quant_refgenome <- function(quant_files, gencode_hla_tx) {

  readDT <- function(x)
    fread(x, select = c(1, 4, 5)
        )[, target_id := sub("^([^|]+).*$", "\\1", target_id)
	][target_id %in% gencode_hla_tx$target_id]

  ldt <- parallel::mclapply(quant_files, readDT, mc.cores = 50)

  dt <- 
    rbindlist(ldt, idcol = "subject")[, .(subject, target_id, est_counts, tpm)]

  dt[gencode_hla_tx, on = .(target_id), nomatch = 0L
   ][, .(subject, gene_id, gene_name, target_id, est_counts, tpm)
   ][order(subject, gene_name)]
}

samples <- sprintf("sample_%02d", 1:50)

quant_round <- commandArgs(TRUE)[1]
quant_dir <- paste0("./quantifications_", quant_round)
quant_files <- file.path(quant_dir, samples, "abundance.tsv")
names(quant_files) <- samples
  
if (quant_round == 1 | quant_round == 2) {
  
  imgt_quants <- process_quant_imgt(quant_files)
  
  hla_names <- imgt_quants[, .(target_id = unique(target_id))]
  hla_names[, `:=`(locus = sub("HLA-", "", imgt_to_gname(target_id)),
		   gene_id = imgt_to_gid(target_id))]

  imgt_quants <- 
    imgt_quants[hla_names, on = .(target_id)
	      ][, .(subject, gene_id, locus, allele = target_id, est_counts, tpm)]
  
  if (quant_round == 1) {
     
    true_genos <- 
      fread("../../data/genos.tsv")[, locus := sub("HLA-", "", locus)]
 
    thresholds <- as.list(seq(0, .25, .05))
    names(thresholds) <- seq(0, .25, .05)

    typings <-
      parallel::mclapply(thresholds, 
			 function(th) hla_genotype_dt(imgt_quants, th),
			 mc.cores = length(thresholds))
      
    typings_dt <- rbindlist(typings, idcol = "th")

    genos_dt <- typings_dt[locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1"), 
			   .(th, subject, locus, allele)]
    genos_dt[, allele := hla_trimnames(gsub("IMGT_|_s\\d", "", allele), 3)]

    accuracies_list <- 
      parallel::mclapply(split(genos_dt, genos_dt$th),
			 function(i) calc_genotyping_accuracy(i, true_genos),
			 mc.cores = length(thresholds))

    accuracies_dt <- rbindlist(accuracies_list, idcol = "th")
    accuracies_dt[, th_average := mean(accuracy), by = th] 

    best_th <- accuracies_dt[which.max(th_average), th]

    accuracies_dt[, `:=`(accuracy = round(accuracy, 2), 
			 th_average = round(th_average, 3))]

    fwrite(accuracies_dt, "./genotyping_accuracies.tsv", 
	   sep = "\t", quote = FALSE)

    out_dt <- typings_dt[th == best_th]
    out_dt[, th := NULL]

  } else if (quant_round == 2) {

  out_dt <- hla_genotype_dt(imgt_quants, th = 0)

  }

} else if (quant_round == "CHR" | quant_round == "ALL") {

  hla_genes <- paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))

  gencode_hla_tx <- 
    setDT(gencode_all_tx
	)[gene_name %in% hla_genes, .(target_id = tx_id, gene_id, gene_name)]

  by_tx <- process_quant_refgenome(quant_files, gencode_hla_tx)

  by_gene <- by_tx[, .(est_counts = sum(est_counts), tpm = sum(tpm)), 
		   by = .(subject, gene_name)]
	
  out_dt <- 
    by_gene[, locus := sub("HLA-", "", gene_name)
	  ][order(subject, locus), .(subject, locus, est_counts, tpm)]

} else {

  stop("wrong dir")
 
}

fwrite(out_dt, file.path(quant_dir, "processed_quant.tsv"), 
       sep = "\t", quote = FALSE)
