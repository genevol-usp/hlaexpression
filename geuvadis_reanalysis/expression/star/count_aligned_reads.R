library(tidyverse)

doMC::registerDoMC(50)

samples <- readLines("../../data/sample_info/samples_phase3_ena_eur.txt") 

aligned_reads <-
  paste0("./quantifications_2/", samples, "/logs/salmon_quant.log") %>%
  setNames(samples) %>%
  plyr::ldply(. %>% 
	      readLines() %>%
	      grep("^Total # of mapped reads", ., value = TRUE) %>%
	      sub("Total # of mapped reads : (\\d+)$", "\\1", .) %>%
	      as.integer(),
	    .id = "subject", .parallel = TRUE)

write_tsv(aligned_reads, "./aligned_reads.tsv")
