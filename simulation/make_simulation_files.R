library(data.table)
library(Biostrings)
devtools::load_all("~/genomicRutils")

### BED

geuvadis_dt <- 
  setDT(geuvadis_info)[kgp_phase3 == 1, .(sample_id = name, subject = ena_id)]

abundance_files <- 
  file.path("../quantifications_2", geuvadis_dt$subject, "abundance.tsv")
names(abundance_files) <- geuvadis_dt$subject

abundances <- 
  parallel::mclapply(abundance_files, function(i) fread(i, select = c(1, 4, 5)),
		     mc.cores = 50)

abundances_dt <- 
  rbindlist(abundances, idcol = "subject")[
  geuvadis_dt, on = "subject"]

abundances_wide <- 
  dcast(abundances_dt, target_id ~ sample_id, value.var = "est_counts")
 
samples <- names(abundances_wide)[-1]

abundances_wide[, (samples) := lapply(.SD, function(x) round(ifelse(is.na(x), 0, x))),
		.SDcols = samples]

fwrite(abundances_wide, "./ground_truth_files/phenotypes.tsv", 
       sep = "\t", quote = FALSE)
