samplesDT <- data.table::fread("../data/nci_sample_info.tsv")

writeLines(unique(samplesDT$sample_id), "./nci_samples.txt")
