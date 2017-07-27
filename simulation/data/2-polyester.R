library(polyester)

counts_matrix <- data.matrix(data.table::fread("./phenotypes.tsv", drop = 1))

simulate_experiment_countmat(fasta = "./polyester_index.fa", 
			     readmat = counts_matrix,
			     readlen = 75, 
			     outdir = "./fastq")
