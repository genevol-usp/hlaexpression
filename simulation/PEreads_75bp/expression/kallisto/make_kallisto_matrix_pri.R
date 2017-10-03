devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

gencode <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id)]

samples <- sprintf("sample_%02d", 1:50) 

abundance_files <- file.path("./quantifications_PRI", samples, "abundance.tsv")
names(abundance_files) <- samples

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject"
	  )[, .(subject, target_id, tpm)
	  ][gencode, on = .(target_id), nomatch = 0L
	  ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- 
  expression_dt[, .(mean(gene_tpm > 0)), by = .(gene_id)][V1 >= 0.9]

expression_dt <- expression_dt[gene_id %in% expressedGenes$gene_id]

expression_w <- dcast(expression_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(expression_w, "./kallisto_PRI_expressed90%.csv")
