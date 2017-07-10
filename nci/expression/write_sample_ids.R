
samplesDT <- data.table::fread("./nci_sample_info.tsv")

writeLines(unique(samplesDT$sample_id), "./nci_samples.txt")
