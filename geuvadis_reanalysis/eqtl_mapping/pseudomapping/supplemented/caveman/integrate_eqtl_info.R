devtools::load_all("~/hlaseqlib")
library(tidyverse)

qtls <- read_qtltools("./conditional_eqtl_best_significant.txt") %>%
    select(phen_id, var_chr, var_id) %>%
    mutate(var_chr = as.integer(var_chr))

info <- read_delim("./eqtl_info.txt", col_names = FALSE, delim = " ")

info %>% 
    left_join(qtls, by = c("X1" = "var_chr", "X3" = "var_id")) %>%
    select(phen_id, X1, X2, X4, X5) %>%
    write_tsv("./eqtl.list", col_names = FALSE)

bed <- read_tsv("./phenotypes.bed")

bed %>%
    filter(id %in% qtls$phen_id) %>%
    write_tsv("./phenotypes_with_eqtl.bed")
