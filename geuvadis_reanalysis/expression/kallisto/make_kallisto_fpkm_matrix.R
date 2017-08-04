devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(geuvadis_info)
setDT(gencode_chr_tx)
setDT(gencode_chr_gene)

gencode <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id, gene_name)]

samples_dt <- geuvadis_info[kgp_phase3 == 1, .(name, subject = ena_id)]

quant_files <-
  file.path("./quantifications_2", samples_dt$subject, "abundance.tsv")
names(quant_files) <- samples_dt$subject

abundances_list <- parallel::mclapply(quant_files, fread, mc.cores = 50)

abundances_dt <- 
  rbindlist(abundances_list, idcol = "subject"
	  )[samples_dt, on = .(subject)
	  ][, fpkm := counts_to_fpkm(est_counts, eff_length), by = subject
          ][, .(subject = name, target_id, fpkm)]

target_set <- abundances_dt[, .(target_id = unique(target_id))]

autosomes_set <- 
  target_set[gencode, on = .(target_id), nomatch = 0L][, .(target_id, gene_id)]

imgt_set <-
  target_set[grepl("^IMGT", target_id)
	   ][, gene_name := imgt_to_gname(target_id)
	   ][gencode, on = .(gene_name), nomatch = 0L
	   ][order(target_id), .(target_id, gene_id)]

gene_set <- rbind(autosomes_set, imgt_set)

gene_dt <- abundances_dt[gene_set, on = .(target_id), nomatch = 0L
		       ][, .(gene_fpkm = sum(fpkm)), by = .(subject, gene_id)]

expressed_genes <- gene_dt[, .(mean(gene_fpkm > 0)), by = .(gene_id)][V1 >= 0.9]

gene_dt <- gene_dt[gene_id %in% expressed_genes$gene_id]

gene_dt_wide <- dcast(gene_dt, subject ~ gene_id, value.var = "gene_fpkm")

fwrite(gene_dt_wide, "./kallisto_gene_expressed90%_fpkm.csv")
