library(data.table)
devtools::load_all("~/genomicRutils")

samples <- geuvadis_info$ena_id[geuvadis_info$kgp_phase3 == 1]

quant_files <-
  file.path("./quantifications_2", samples, "abundance.tsv")
names(quant_files) <- samples

abundances <- 
  rbindlist(parallel::mclapply(quant_files, fread, mc.cores = 50), 
	    idcol = "subject")

abundances[, fpkm := counts_to_fpkm(est_counts, eff_length), by = subject]

hla <- abundances[grepl("^IMGT_(A|B|C|DQA1|DQB1|DRB1)\\*", target_id)]

hla[, locus := imgt_to_gname(target_id)]

fpkm_per_hla <- hla[, .(fpkm = sum(fpkm)), by = .(subject, locus)]

fwrite(fpkm_per_hla, "./quantifications_2/gene_fpkms.tsv",
       sep = "\t", quote = FALSE)
