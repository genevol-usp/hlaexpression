library(tidyverse)

pcr <- read_tsv("./hla_carrington_data_filteredN10.tsv") %>%
    mutate(lineage = sub("\\*", "", lineage))

observed_means <- pcr %>%
    group_by(lineage) %>%
    summarise(m = mean(rna)) %>%
    ungroup()

lm_df <- pcr %>%
    select(-rna) %>%
    mutate(hit = 1) %>%
    complete(subject, nesting(locus, lineage), fill = list(hit = 0)) %>%
    select(subject, locus, lineage, hit) %>%
    arrange(subject, lineage) %>%
    group_by(subject, locus, lineage) %>%
    summarise(hit = sum(hit)) %>%
    ungroup() %>%
    inner_join(distinct(pcr, subject, locus, rna), by = c("subject", "locus"))

estimate_df <- lm_df %>%
    split(.$locus) %>%
    map(~spread(., lineage, hit)) %>%
    map(~select(., -subject, -locus)) %>%
    lapply(function(x) broom::tidy(lm(rna ~ . - 1, data = x))) %>%
    map(~select(., 1:2)) %>%
    bind_rows(.id = "locus") %>%
    select(locus, lineage = term, y = estimate)

fitted_df <- lm_df %>%
    split(.$locus) %>%
    map(~spread(., lineage, hit)) %>%
    map(~mutate(., f = paste("rna ~", paste(names(.)[-(1:3)], collapse = "+"), "- 1"),
		y = fitted(lm(as.formula(f))))) %>%
    map(~select(., subject, locus, rna, y)) %>%
    bind_rows() %>%
    left_join(select(pcr, subject, locus, lineage), by = c("subject", "locus")) %>%
    group_by(subject, locus, rna, y) %>%
    summarise(lineage = paste(lineage, collapse = "/")) %>%
    ungroup() %>%
    select(subject, locus, lineage, rna, y)

write_tsv(fitted_df, "./carrington_data_lmEstimates.tsv")



### shrinkage correction

estimate_df %>% group_by(locus) %>% summarise(mean(y))

se <- fitted_df %>% filter(locus == "A") %>% summarise(se = sd(y)/n()) %>% pull(se)

fitted_df %>% filter(locus == "A") %>%
    mutate(adj_y = y / (10000000 - (1 * se))) %>%
    select(subject, adj_y) %>%
    left_join(filter(lm_df, locus == "A"), .) %>%
    select(subject, locus, lineage, hit, rna = adj_y) %>%
    spread(lineage, hit) %>%
    split(.$locus) %>%
    map(~select(., -subject, -locus)) %>%
    lapply(function(x) broom::tidy(lm(rna ~ . - 1, data = x)))




###


means_df <- left_join(observed_means, estimate_df, by = "lineage") %>%
    select(locus, lineage, obs = m, exp = y)

x <- dat_filtered %>%
    left_join(means_df, by = c("locus", "lineage")) %>%
    group_by(subject, locus) %>%
    summarise(rna = unique(rna), obs = sum(obs), exp = sum(exp)) %>%
    ungroup()

summary(lm(obs ~ rna, data = filter(x, locus == "A")))
broom::tidy(lm(obs ~ rna, data = filter(x, locus == "A")))

summary(lm(exp ~ rna, data = filter(x, locus == "A")))
broom::tidy(lm(exp ~ rna, data = filter(x, locus == "A")))

