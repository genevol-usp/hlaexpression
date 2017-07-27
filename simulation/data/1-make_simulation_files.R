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

writeLines(geuvadis_dt$sample_id, "./samples.txt")

abundance_files <- 
  file.path("../../geuvadis_reanalysis/expression/kallisto/quantifications_2", 
	    sample_dt$subject, "abundance.tsv")
names(abundance_files) <- sample_dt$subject

abundances <- 
  parallel::mclapply(abundance_files, function(i) fread(i, select = c(1, 4, 5)),
		     mc.cores = 25)

abundances_dt <- rbindlist(abundances, idcol = "subject")[sample_dt, on = "subject"]

genos <- abundances_dt[grepl("^IMGT_(A|B|C|DQA1|DQB1|DRB1)", target_id) & est_counts > 0]
genos[, `:=`(locus = imgt_to_gname(target_id), 
	     allele = hla_trimnames(gsub("IMGT_", "", target_id), 3))]
hom <- copy(genos)[, n := .N, by = .(subject, locus)][n == 1L][, n := NULL]
genos <- rbind(genos, hom)[, .(subject = code, locus, allele)][order(subject, locus, allele)] 

fwrite(genos, "./genos.tsv", sep = "\t", quote = FALSE)

index <- readDNAStringSet("../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.IMGT.transcripts.fa")
index <- index [width(index) >= 75]

abundances_wide <- dcast(abundances_dt, target_id ~ sample_id, value.var = "est_counts")

tx <- intersect(names(index), abundances_wide$target_id)

index <- index[tx]
abundances_wide <- abundances_wide[target_id %in% tx]

samples <- names(abundances_wide)[-1]

abundances_wide[, (samples) := lapply(.SD, function(x) round(ifelse(is.na(x), 0, x))), .SDcols = samples]
abundances_wide[, (samples) := lapply(.SD, function(x) round(x/sum(x) * 3e7)), .SDcols = samples]

writeXStringSet(index, "./polyester_index.fa")
fwrite(abundances_wide, "./phenotypes.tsv", sep = "\t", quote = FALSE)
