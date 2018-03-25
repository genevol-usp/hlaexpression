library(data.table)

pcs <- fread("./eur.pca", nrows = 3)
pcs[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]

setnames(pcs, "SampleID", "id")

fwrite(pcs, "./covariates_genos.txt", sep = "\t", quote = FALSE)
