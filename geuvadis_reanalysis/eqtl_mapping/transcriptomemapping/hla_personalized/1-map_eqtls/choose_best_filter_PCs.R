library(tidyverse)

pcs <- seq(0, 100, 10)
th <- c(25, 50, 75, 90)

egenes_df <- 
    expand.grid(th, pcs) %>%
    mutate(path = sprintf("./th_%d/2-permutations/results/permutations_%d.significant.txt", Var1, Var2),
	   eGenes = map_int(path, ~nrow(read_delim(., col_names = FALSE, delim = " ")))) %>%
    rename(th = Var1, pc = Var2) %>%
    arrange(desc(eGenes))

write_tsv(egenes_df, "./best_filter_PCs.tsv")
