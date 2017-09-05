library(data.table)
library(Biostrings)
devtools::load_all("/home/vitor/hlaseqlib")
  
setDT(geuvadis_info)

set.seed(2)

sample_dt <- 
  geuvadis_info[kgp_phase3 == 1, .(sample_id = name, subject = ena_id)
	      ][sample_id %in% sample(sample_id, 50)
	      ][order(sample_id)
	      ][, code := sprintf("sample_%02d", 1:50)]

writeLines(sample_dt$sample_id, "./samples.txt")

abundance_files <- 
  file.path("../../geuvadis_reanalysis/expression/star/quantifications_2", 
	    sample_dt$subject, "quant.sf")
names(abundance_files) <- sample_dt$subject

abundances <- 
  parallel::mclapply(abundance_files, function(i) fread(i, select = c(1, 4, 5)),
		     mc.cores = 25)

abundances_dt <- 
  rbindlist(abundances, idcol = "subject")[sample_dt, on = "subject"]

genos <- 
  abundances_dt[grepl("^IMGT_(A|B|C|DQA1|DQB1|DRB1)", Name) & NumReads > 0]

genos[, `:=`(locus = imgt_to_gname(Name), 
	     allele = hla_trimnames(gsub("IMGT_", "", Name), 3))]

hom <- copy(genos)[, n := .N, by = .(subject, locus)][n == 1L][, n := NULL]

genos <- 
  rbind(genos, hom
      )[, .(subject = code, locus, allele)
      ][order(subject, locus, allele)] 

fwrite(genos, "./genos.tsv", sep = "\t", quote = FALSE)

index <- 
  readDNAStringSet("../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

abundances_wide <- dcast(abundances_dt, Name ~ sample_id, value.var = "NumReads")
setkey(abundances_wide, Name)

tx <- intersect(names(index), abundances_wide$Name)

index <- index[tx]
abundances_wide <- abundances_wide[tx]

samples <- names(abundances_wide)[-1]

abundances_wide[, (samples) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = samples]
abundances_wide[, (samples) := lapply(.SD, function(x) round(x/sum(x) * 3e7)), .SDcols = samples]

writeXStringSet(index, "./polyester_index.fa")
fwrite(abundances_wide[, -1], "./phenotypes.tsv", sep = "\t", quote = FALSE)
