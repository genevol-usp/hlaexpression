library(tidyverse)

# can't use read table functions because file has different number of columns
regdb <- readLines("./RegulomeDB.dbSNP141.chr6.txt") %>%
    strsplit("\t")

db_df <- tibble(var = map_chr(regdb, 4),
		score = map_chr(regdb, last))

rsmerge <- 
    read_tsv("../../../../../data/previous_qtls/RsMergeArch.bcp.gz", 
	     col_names = FALSE) %>%
    select(rsHigh = X1, rsLow = X2, rsCurrent = X7) %>%
    mutate_at(vars(rsHigh, rsLow, rsCurrent), function(x) paste("rs", x))

db_eqtl_df <- filter(db_df, score %in% paste0(1, letters[1:6])) %>%
    left_join(rsmerge, by = c("var" = "rsHigh")) %>%
    mutate(var = ifelse(is.na(rsCurrent), var, rsCurrent)) %>%
    select(var, score)

write_tsv(db_eqtl_df, "./variants.tsv")
