library(tidyverse)

doMC::registerDoMC(50)

samples <- readLines("../../data/sample_info/samples_phase3_ena_eur.txt") 

aligned_reads <-
  paste0("./quantifications_2/log/", samples, ".quant.log") %>%
  setNames(samples) %>%
  plyr::ldply(. %>% 
	      readLines() %>%
	      grep("reads pseudoaligned$", ., value = TRUE) %>%
	      sub("^\\[quant\\] processed [0-9,]+ reads, ([0-9,]+) reads pseudoaligned$", "\\1", .) %>%
	      gsub(",", "", .) %>%
	      as.integer(),
	    .id = "subject", .parallel = TRUE)

write_tsv(aligned_reads, "./aligned_reads.tsv")
