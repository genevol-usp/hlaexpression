library(data.table)

pcs_dt <- 
  fread("./phenotypes_pcs.pca", nrows = 10
      )[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]
setnames(pcs_dt, "SampleID", "id")

fwrite(pcs_dt, "./covariates_pheno_10.txt", sep = "\t", quote = FALSE)
