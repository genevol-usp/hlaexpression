devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

mismatch <- commandArgs(TRUE)[1]

read_quant <- function(x) {

  fread(x)[, Name := sub("^([^|]+).*$", "\\1", Name)]
}

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

autosomes <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id)]

samples <- sprintf("sample_%02d", 1:50) 

abundance_files <- 
  file.path(paste0("./quantifications_CHR_", mismatch), samples, "quant.sf")
names(abundance_files) <- samples

expression_list <- parallel::mclapply(abundance_files, read_quant, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject"
	  )[, .(subject, target_id = Name, tpm = TPM)
          ][autosomes, on = .(target_id), nomatch = 0L
	  ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- 
  expression_dt[, .(mean(gene_tpm > 0)), by = .(gene_id)][V1 >= 0.9]

expression_dt <- expression_dt[gene_id %in% expressedGenes$gene_id]

out <- dcast(expression_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(out, paste0("./quantifications_CHR_", mismatch, "_expressed90%.csv"))

