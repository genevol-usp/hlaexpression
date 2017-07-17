library(data.table)

pcs <- fread("./pca/eur.pca", nrows = 3)
pcs[, SampleID := sub("^.+(PC\\d+)$", "\\1", SampleID)]

setnames(pcs, "SampleID", "id")

fwrite(pcs, "./pca_genotypes/covariates_genos.txt", sep = "\t", quote = FALSE)
