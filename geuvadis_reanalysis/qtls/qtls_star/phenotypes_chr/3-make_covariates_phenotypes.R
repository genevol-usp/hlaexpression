library(data.table)
devtools::load_all("/home/vitor/hlaseqlib")

pcs_dt <- 
  fread("./phenotypes_eur_pcs.pca", nrows = 10
      )[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]
setnames(pcs_dt, "SampleID", "id")

cols <- names(pcs_dt)[-1]
pcs_dt[, (cols) := lapply(.SD, as.character), .SDcols = cols]

fwrite(pcs_dt, "./covariates/covariates_pheno_10.txt", sep = "\t", quote = FALSE)
