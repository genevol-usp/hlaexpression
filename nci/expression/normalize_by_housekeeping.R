library(data.table)
devtools::load_all("~/genomicRutils")

readDT <- function(i) {
  x <- fread(i, select = c(1, 4, 5))
  x[target_id %in% housekeep$target_id]
}

samples <- readLines("./nci_samples.txt")
quant_files <- file.path("./quantifications_2", samples, "abundance.tsv")
names(quant_files) <- samples

housekeep <- 
  setDT(gencode_chr_tx
      )[gene_name %in% c("B2M", "ACTB", "GAPDH")
      ][, .(gene_name, target_id = tx_id)]

quants <- 
  rbindlist(parallel::mclapply(quant_files, readDT, mc.cores = 48), idcol = "subject")

quants_hk <- 
  quants[housekeep, on = .(target_id)
       ][, .(est_counts = sum(est_counts)), by = .(subject, gene_name)] 

quants_hk <- dcast(quants_hk, subject ~ gene_name, value.var = "est_counts")

quants_imgt <- 
  fread("./quantifications_2/processed_quant.tsv"
      )[locus %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")
      ][quants_hk, on = .(subject)
      ][, `:=`(est_counts_ACTBnorm = est_counts/ACTB,
	       est_counts_B2Mnorm = est_counts/B2M,
	       est_counts_GAPDHnorm = est_counts/GAPDH)
       ][, .(subject, locus, allele, tpm, 
	     est_counts_ACTBnorm, est_counts_B2Mnorm, est_counts_GAPDHnorm)]

fwrite(quants_imgt, "./quantifications_2/housekeeping_norm_expression.tsv",
       sep = "\t", quote = FALSE)
