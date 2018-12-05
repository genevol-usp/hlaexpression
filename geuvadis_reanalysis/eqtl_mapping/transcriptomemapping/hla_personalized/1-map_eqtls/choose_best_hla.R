library(tidyverse)

df <- 
    tibble(th = c(25L, 50L, 75L, 90L),
	   path = sprintf("./th_%i/3-conditional_analysis/hla_qtls.tsv", th),
	   d = map(path, ~filter(read_tsv(.), best == 1))) %>%
    select(-path) %>%
    unnest(d)

df %>% filter(th == 50)

