library(data.table)

perm_files <- 
  sprintf("./results/permutations_%d.significant.txt", seq(0, 100, 5))
names(perm_files) <- seq(0, 100, 5)  
 
perm_dt <- rbindlist(lapply(perm_files, fread), idcol = "PCs")

n_genes <- perm_dt[, .N, by = PCs]

best <- n_genes[order(-N)]
