library(polyester)

icol <- as.integer(commandArgs(TRUE)[1])

counts_matrix <- 
  data.matrix(data.table::fread("./phenotypes.tsv", select = icol))

simulate_experiment_countmat(fasta = "./polyester_index.fa", 
			     readmat = counts_matrix,
			     readlen = 75, 
			     outdir = paste0("./fasta_", icol))
