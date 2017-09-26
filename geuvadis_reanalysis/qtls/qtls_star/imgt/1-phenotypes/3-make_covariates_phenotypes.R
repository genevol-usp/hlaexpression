library(data.table)

pcs_dt <- 
    fread("./phenotypes_eur_pcs.pca", nrows = 100
	)[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]
setnames(pcs_dt, "SampleID", "id")

out_basename <- "./covariates/covariates_pheno_"

for (pc in c(seq(5, 30, 5), seq(40, 100, 10))) {

    dti <- pcs_dt[id %in% paste0("PC", seq_len(pc))]
    fwrite(dti, paste0(out_basename, pc, ".txt"), sep = "\t", quote = FALSE)
}
