library(data.table)
library(Biostrings)
devtools::load_all("/home/vitor/hlaseqlib")
  
setDT(geuvadis_info)

set.seed(2)
sample_dt <- 
  geuvadis_info[kgp_phase3 == 1 & pop != "YRI", .(sample_id = name, subject = ena_id)
	      ][sample_id %in% sample(sample_id, 50)
	      ][order(sample_id)
	      ][, code := sprintf("sample_%02d", 1:50)]

abundance_files <- 
  file.path("../../geuvadis_reanalysis/expression/star/quantifications_2", 
	    sample_dt$subject, "quant.sf")
names(abundance_files) <- sample_dt$subject

abundances <- 
  parallel::mclapply(abundance_files, function(i) fread(i, select = c(1, 4, 5)),
		     mc.cores = 25)

abundances_imgt <- 
    fread("../../geuvadis_reanalysis/expression/star/quantifications_2/processed_quant.tsv"
	)[sample_dt, on = .(subject)
	][, .(subject = code, locus, Name = allele, NumReads = est_counts)]

genos <- 
    abundances_imgt[locus %in% paste0("HLA-", c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1"))
		  ][, .(subject, locus, allele = Name)]
genos[, allele := hla_trimnames(gsub("IMGT_", "", sub("^([^/]+).*$", "\\1", allele)), 3)] 

fwrite(genos, "./genos.tsv", sep = "\t", quote = FALSE)

abundances_imgt[, locus := NULL]

abundances_non_imgt <- 
  rbindlist(abundances, idcol = "subject"
	  )[sample_dt, on = "subject"
	  ][, .(subject = code, Name, NumReads)
	  ][!grepl("^IMGT", Name)]

abundances_dt <-
    rbind(abundances_non_imgt, abundances_imgt
	)[, .(NumReads = sum(NumReads)), by = .(subject, Name)]

index <- readDNAStringSet("../../imgt_index/gencode.v25.PRI.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

abundances_wide <- dcast(abundances_dt, Name ~ subject, value.var = "NumReads")
setkey(abundances_wide, Name)

tx <- intersect(names(index), abundances_wide$Name)

index <- index[tx]
abundances_wide <- abundances_wide[tx]

samples <- sample_dt$code

abundances_wide[, (samples) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = samples]
abundances_wide[, (samples) := lapply(.SD, function(x) round(x/sum(x) * 3e7)), .SDcols = samples]

writeXStringSet(index, "./polyester_index.fa")
fwrite(abundances_wide[, -1], "./phenotypes.tsv", sep = "\t", quote = FALSE)
