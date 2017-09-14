devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

gencode <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id, gene_name)]

samples <- sprintf("sample_%02d", 1:50)

abundance_files <- file.path("./quantifications_2", samples, "quant.sf")
names(abundance_files) <- samples

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject"
	  )[, .(subject, target_id = Name, tpm = TPM)]

target_set <- expression_dt[, .(target_id = unique(target_id))]

autosomes_set <- 
  target_set[gencode, on = .(target_id), nomatch = 0L][, .(target_id, gene_id)]

imgt_set <-
  target_set[grepl("^IMGT", target_id)
	   ][, gene_name := imgt_to_gname(target_id)
	   ][gencode, on = .(gene_name), nomatch = 0L
	   ][order(target_id), .(target_id, gene_id)]

gene_set <- rbind(autosomes_set, imgt_set)

gene_dt <- expression_dt[gene_set, on = .(target_id), nomatch = 0L
		       ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- gene_dt[, .(mean(gene_tpm > 0)), by = .(gene_id)][V1 >= 0.9]

gene_dt <- gene_dt[gene_id %in% expressedGenes$gene_id]

gene_dt_wide <- dcast(gene_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(gene_dt_wide, "./quantifications_imgt_expressed90%.csv")
