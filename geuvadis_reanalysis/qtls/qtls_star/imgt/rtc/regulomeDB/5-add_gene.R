library(tidyverse)

extract_genes <- function(x) {

    x_split <- strsplit(x, ",") %>% 
	unlist() 
    
    eQTL <- x_split %>% 
	grep("eQTL", ., value = TRUE) %>%
	strsplit("\\|") %>% 
	map(3) %>% 
	unlist()

    dsQTL <- x_split %>%
	grep("dsQTL", ., value = TRUE) %>%
	strsplit("\\|") %>% 
	map(5) %>% 
	unlist()

    paste(c(eQTL, dsQTL), collapse = ",")
}

# can't use read table functions because file has different number of columns
# can't fill missing columns either, bc they are in different order across rows!
regdb <- readLines("./RegulomeDB.dbSNP141.chr6.txt") %>%
    strsplit("\t")

db_df <- tibble(var = map_chr(regdb, 4),
		score = map_chr(regdb, last))

results <- read_tsv("./results.tsv")

# works for this subset of variants, but other subsets may need to update rs ids
# or the list elements may have different lengths
reg_subset <- regdb[which(db_df$var %in% results$rtc_var)]

affected_genes <- reg_subset %>% map(26) %>% map_chr(extract_genes)

genes_df <- tibble(rtc_var = map_chr(reg_subset, 4),
		   affected_gene = affected_genes)

left_join(results, genes_df) %>% 
    rename(score = info) %>%
    group_by(gene, variant, rank, d_prime, rtc, score, affected_gene) %>%
    summarize(rtc_var = paste(rtc_var, collapse = ",")) %>%
    ungroup() %>%
    select(gene, rank, variant, rtc_var, d_prime, rtc, score, affected_gene) %>%
    mutate(i = map2_lgl(gene, affected_gene, function(x, y) grepl(x, y, fixed = TRUE))) %>%
    group_by(gene, rank, variant, i) %>%
    filter(rtc == max(rtc)) %>%
    ungroup() %>% 
    arrange(gene, rank, desc(rtc)) %>%
    select(-i) %>%
    write_tsv("./results_affected_genes.tsv")
