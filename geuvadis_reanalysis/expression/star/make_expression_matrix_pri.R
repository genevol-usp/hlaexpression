devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

autosomes <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id)]

samples_dt <- 
  geuvadis_info[kgp_phase3 == 1 & pop != "YRI", .(name, subject = ena_id)]

abundance_files <- 
    file.path("./quantifications_PRI", samples_dt$subject, "quant.sf")
names(abundance_files) <- samples_dt$subject

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject"
	  )[samples_dt, on = .(subject)
	  ][, .(subject = name, target_id = Name, tpm = TPM)
          ][autosomes, on = .(target_id), nomatch = 0L
	  ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

expressedGenes <- 
  expression_dt[, .(mean(gene_tpm >= 1)), by = .(gene_id)][V1 >= 0.9]

expression_dt <- expression_dt[gene_id %in% expressedGenes$gene_id]

out <- dcast(expression_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(out, "./quantifications_pri_expressed90%.csv")

