library(data.table)

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

perm_files <- sprintf("./results/permutations_%d.significant.txt", pcs)
names(perm_files) <- pcs
 
perm_dt <- rbindlist(lapply(perm_files, fread), idcol = "PCs")

n_genes <- perm_dt[, .N, by = PCs]

n_genes[order(-N)]
