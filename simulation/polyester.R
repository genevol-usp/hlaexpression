library(data.table)
library(polyester)
library(Biostrings)
devtools::load_all("~/hlatools")
devtools::load_all("~/genomicRutils")

index <- readDNAStringSet("../index/gencode.v26.CHR.IMGT.transcripts.fa")
index <- index[width(index) >= 75]

ground_truth <- fread("./ground_truth_files/phenotypes.tsv")

set.seed(2)
all_samples <- names(ground_truth)[-1]
chosen_samples <- sort(sample(all_samples, 50))

samples_dt <- data.table(subject = chosen_samples,
			 code = sprintf("sample_%02d", 1:50))

writeLines(chosen_samples, "./ground_truth_files/samples.txt")

counts_dt <- ground_truth[, c("target_id", (chosen_samples)), with = FALSE]
setkey(counts_dt, target_id)

genos <- 
  melt(counts_dt[grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", target_id)], 
       id = 1, measure = 2:ncol(counts_dt), variable = "subject"
      )[value > 0
      ][, locus := imgt_to_gname(target_id)
      ][, allele := hla_trimnames(gsub("IMGT_|_s\\d", "", target_id), 3)
      ][, .(subject, locus, allele)]

hom <- copy(genos)[, n := .N, by = .(subject, locus)][n == 1][, n := NULL]

genos <- 
  rbind(genos, hom)[order(subject, locus, allele)
		  ][samples_dt, on = .(subject)
		  ][, .(subject = code, locus, allele)]

fwrite(genos, "./ground_truth_files/genos.tsv", sep = "\t", quote = FALSE)

tx <- intersect(names(index), counts_dt$target_id)

index <- index[tx]

counts_dt <- 
  counts_dt[tx][
  , lapply(.SD, function(x) round(x/sum(x) * 3e7)), .SDcols = chosen_samples]

counts_matrix <- as.matrix(counts_dt)

writeXStringSet(index, "./ground_truth_files/polyester_index.fa")
  
simulate_experiment_countmat(fasta = "./ground_truth_files/polyester_index.fa", 
			     readmat = counts_matrix,
			     readlen = 75, 
			     outdir = "./fastq")
