devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

setDTthreads(36)
setDT(gencode_chr_tx)

quant_round <- commandArgs(TRUE)[1]
quant_dir <- paste0("./quantifications_", quant_round)

gencode <- gencode_chr_tx[chr %in% 1:22, .(target_id = tx_id, gene_id, gene_name)]

samples <- sprintf("sample_%02d", 1:50)

samples_dt <- data.table(subject = readLines("../../data/samples.txt"),
			 code = samples)

index <- Biostrings::readDNAStringSet("../../data/polyester_index.fa")

imgt_index_ids <- data.table(target_id = grep("IMGT", names(index), value = TRUE))
imgt_index_ids[, gene_id := imgt_to_gid(target_id)]

all_index_ids <- 
  rbind(gencode_chr_tx[, .(target_id = tx_id, gene_id)], imgt_index_ids)

ground_exp <- fread("../../data/phenotypes.tsv")
ground_exp <- 
  ground_exp[, target_id := names(index)
	   ][all_index_ids, on = .(target_id), nomatch = 0L]

ground_exp <- 
  melt(ground_exp, id = c("target_id", "gene_id"), measure.vars = 1:50,
       variable.name = "subject", value.name = "counts")

ground_exp <- 
  ground_exp[, .(counts = sum(counts)), by = .(subject, gene_id)
	   ][samples_dt, on = .(subject)][, .(subject = code, gene_id, counts)
	   ][!is.na(gene_id)]

abundance_files <- file.path(quant_dir, samples, "abundance.tsv")
names(abundance_files) <- samples

expression_list <- parallel::mclapply(abundance_files, fread, mc.cores = 50)

expression_dt <- 
  rbindlist(expression_list, idcol = "subject")[, .(subject, target_id, est_counts)]

target_set <- expression_dt[, .(target_id = unique(target_id))]

autosomes_set <- 
  target_set[gencode, on = .(target_id), nomatch = 0L][, .(target_id, gene_id)]

imgt_set <-
  target_set[grepl("^IMGT", target_id)
	   ][, gene_name := imgt_to_gname(target_id)
	   ][gencode, on = .(gene_name), nomatch = 0L
	   ][order(target_id), .(target_id, gene_id)]

gene_set <- rbind(autosomes_set, imgt_set)

gene_dt <- 
  expression_dt[gene_set, on = .(target_id), nomatch = 0L
	      ][, .(gene_counts = sum(est_counts)), by = .(subject, gene_id)]

mrd_dt <- 
  ground_exp[gene_dt, on = .(subject, gene_id), nomatch = 0L
	   ][, rd := ifelse(counts == 0 & gene_counts == 0, 0, 
			    abs(counts - gene_counts)/(counts + gene_counts))
	   ][, .(mrd = mean(rd)), by = .(subject)]

fwrite(mrd_dt, paste0("./mrd_", quant_round, ".tsv"), sep = "\t", quote = FALSE)
