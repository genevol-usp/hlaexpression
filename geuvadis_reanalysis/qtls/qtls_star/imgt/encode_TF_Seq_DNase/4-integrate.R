library(tidyverse)

left_join(read_tsv("./hla.TF.intersect.bed", col_names = FALSE) %>% distinct(),
	  read_tsv("./TF.hla.intersect.bed", col_names = FALSE),
	  by = c("X1", "X2", "X3")) %>%
distinct() %>%
write_tsv("./hla.TF.integrated.bed", col_names = FALSE)

left_join(read_tsv("./hla.DNase.intersect.bed", col_names = FALSE) %>% distinct(),
	  read_tsv("./DNase.hla.intersect.bed", col_names = FALSE),
	  by = c("X1", "X2", "X3")) %>%
distinct() %>%
write_tsv("./hla.DNase.integrated.bed", col_names = FALSE)

left_join(read_tsv("./hla.Seg.intersect.bed", col_names = FALSE) %>% distinct(),
	  read_tsv("./Seg.hla.intersect.bed", col_names = FALSE),
	  by = c("X1", "X2", "X3")) %>%
distinct() %>%
write_tsv("./hla.Seg.integrated.bed", col_names = FALSE)

