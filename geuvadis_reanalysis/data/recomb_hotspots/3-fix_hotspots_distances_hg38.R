library(tidyverse)

hg38 <- 
  read_tsv("./hotspots_hg38.bed", col_names = FALSE) %>%
  filter(X1 %in% paste0("chr", 1:22)) %>%
  mutate(X1 = as.integer(gsub("^chr", "", X1)),
         X4 = abs(X3 - X2)) %>%
  arrange(X1, X2)

write_tsv(hg38, "./hotspots_hg38_fixed_distances.bed", col_names = FALSE)
