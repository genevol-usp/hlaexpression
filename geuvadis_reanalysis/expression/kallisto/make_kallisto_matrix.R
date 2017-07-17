devtools::load_all("~/genomicRutils")
library(data.table)

setDTthreads(36)

samples_dt <- setDT(geuvadis_info)[kgp_phase3 == 1, .(name, subject = ena_id)]

abundance_files <- 
  file.path("./quantifications_2", samples_dt$subject, "abundance.tsv")
names(abundance_files) <- samples_dt$subject

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject"
	  )[samples_dt, on = .(subject)
	  ][, .(subject = name, target_id, tpm)]

gencode_autosomes <- 
  setDT(gencode_chr_tx)[chr %in% 1:22, .(target_id = tx_id, gene_id)]

imgt_ids <- 
  expression_dt[grepl("^IMGT", target_id)
	      ][, .(target_id = unique(target_id))
	      ][, gene_id := imgt_to_gid(target_id)
	      ][order(target_id)]

gene_set <- rbind(gencode_autosomes, imgt_ids)[!is.na(gene_id)]

gene_dt <- 
  expression_dt[gene_set, on = .(target_id), nomatch = 0L
	      ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- gene_dt[, .(mean(gene_tpm > 0)), by = .(gene_id)][V1 >= 0.9]

gene_dt <- gene_dt[gene_id %in% expressedGenes$gene_id]

gene_dt_wide <- dcast(gene_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(gene_dt_wide, "./kallisto_gene_expressed90%.csv")
