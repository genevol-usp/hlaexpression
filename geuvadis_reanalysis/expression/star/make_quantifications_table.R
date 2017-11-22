devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

quant_dir <- commandArgs(TRUE)[1]

gencode <- 
    gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id, gene_name)]

samples_dt <- 
    geuvadis_info[kgp_phase3 == 1 & pop != "YRI", .(name, subject = ena_id)]

abundance_files <- file.path(quant_dir, samples_dt$subject, "quant.sf")
names(abundance_files) <- samples_dt$subject

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)
    
expression_dt <- 
    rbindlist(expression_list, idcol = "subject"
	    )[samples_dt, on = .(subject)
	    ][, .(subject = name, target_id = Name, tpm = TPM)]

if (grepl("quantifications_2", quant_dir)) {

    target_set <- expression_dt[, .(target_id = unique(target_id))]

    autosomes_set <- 
	target_set[gencode, on = .(target_id), nomatch = 0L
		 ][, .(target_id, gene_id)]

    imgt_set <-
	target_set[grepl("^IMGT", target_id)
		 ][, gene_name := imgt_to_gname(target_id)
		 ][gencode, on = .(gene_name), nomatch = 0L
		 ][order(target_id), .(target_id, gene_id)]

    gene_set <- rbind(autosomes_set, imgt_set)

    gene_dt <- 
	expression_dt[gene_set, on = .(target_id), nomatch = 0L
		    ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]

} else if (grepl("quantifications_PRI", quant_dir)) {

    gene_dt <- 
	expression_dt[gencode, on = .(target_id), nomatch = 0L
		    ][, .(gene_tpm = sum(tpm)), by = .(subject, gene_id)]
}

expressedGenes <- gene_dt[, .(mean(gene_tpm >= 1)), by = .(gene_id)][V1 >= 0.9]

gene_dt <- gene_dt[gene_id %in% expressedGenes$gene_id]

gene_dt_wide <- dcast(gene_dt, subject ~ gene_id, value.var = "gene_tpm")

fwrite(gene_dt_wide, paste0(quant_dir, "_expressed90%.csv"))
