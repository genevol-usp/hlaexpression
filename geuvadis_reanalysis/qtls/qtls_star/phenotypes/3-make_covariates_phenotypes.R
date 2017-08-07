library(data.table)
devtools::load_all("/home/vitor/hlaseqlib")

known_covs <- 
  setDT(geuvadis_info
      )[pop != "YRI" & kgp_phase3 == 1
      ][, .(subject = name, lab = lab_code, pop)
      ][, `:=`(pop = substring(pop, 1, 1), lab = LETTERS[lab])]

known_covs <- 
  melt(known_covs, id = 1, measure = c("lab", "pop"), variable.name = "id")

known_covs <- dcast(known_covs, id ~ subject, value = "value")

pcs_dt <- 
  fread("./phenotypes_eur_pcs.pca", nrows = 100)[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]
setnames(pcs_dt, "SampleID", "id")

cols <- names(pcs_dt)[-1]
pcs_dt[, (cols) := lapply(.SD, as.character), .SDcols = cols]

out_basename <- "./covariates/covariates_pheno_"

for (pc in c(seq(0, 30, 5), seq(40, 100, 10))) {

  if (pc == 0) {
    fwrite(known_covs, paste0(out_basename, pc, ".txt"), sep = "\t", quote = FALSE)
  } else {
    dti <- pcs_dt[id %in% paste0("PC", seq_len(pc))]
    outi <- rbind(known_covs, dti)
    fwrite(outi, paste0(out_basename, pc, ".txt"), sep = "\t", quote = FALSE)
  }
}
