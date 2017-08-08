devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

make_table <- function(quant_dir) {
  abundance_files <- file.path(quant_dir, samples_dt$subject, "abundance.tsv")
  names(abundance_files) <- samples_dt$subject

  expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

  expression_dt <- 
    rbindlist(expression_list, idcol = "subject"
	    )[samples_dt, on = .(subject)
	    ][, .(subject = name, target_id, tpm)
	    ][gencode, on = .(target_id), nomatch = 0L
	    ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

  expressedGenes <- 
    expression_dt[, .(mean(gene_tpm > 0)), by = .(gene_id)][V1 >= 0.9]

  expression_dt <- expression_dt[gene_id %in% expressedGenes$gene_id]

  dcast(expression_dt, subject ~ gene_id, value.var = "gene_tpm")
}

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

gencode <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id)]

samples_dt <- 
  geuvadis_info[kgp_phase3 == 1 & pop != "YRI", .(name, subject = ena_id)]

make_table("./quantifications_CHR") %>%
  fwrite("./kallisto_CHR_expressed90%.csv")

make_table("./quantifications_ALL") %>%
  fwrite("./kallisto_ALL_expressed90%.csv")

