library(tidyverse)
library(readxl)

dat <- read_excel("./HLA-A, -B, -C expression levels to Vitor.xlsx") %>%
    select(subject = 1, A = 2, B = 3, C = 4, surface_c = 5,
           allele_A1 = 6, allele_A2 = 7, allele_B1 = 8, allele_B2 = 9,
           allele_C1 = 10, allele_C2 = 11)

rna <- select(dat, subject, A:C) %>%
    gather(locus, rna, -1) %>%
    mutate(rna = as.numeric(rna)) %>%
    filter(!is.na(rna)) %>%
    arrange(subject, locus)

genos <- select(dat, subject, allele_A1:allele_C2) %>%
    gather(locus, allele, -1) %>%
    mutate(locus = sub("^allele_([ABC])[1-2]$", "\\1", locus),
	   lineage = sub("^([^:]+).*$", "\\1", allele),
	   lineage = paste0(locus, lineage)) %>%
    select(subject, locus, lineage)

dat_filtered <- inner_join(genos, rna, by = c("subject", "locus")) %>%
    select(subject, locus, lineage, rna) %>%
    group_by(lineage) %>%
    mutate(nLineage = n()) %>%
    group_by(subject, locus) %>%
    filter(all(nLineage >= 5)) %>%
    ungroup() %>%
    select(-nLineage) %>%
    arrange(subject, locus, lineage)

observed_means <- dat_filtered %>%
    group_by(lineage) %>%
    summarise(m = mean(rna)) %>%
    ungroup()

lm_df <- dat_filtered %>%
    select(-rna) %>%
    mutate(hit = 1) %>%
    complete(subject, nesting(locus, lineage), fill = list(hit = 0)) %>%
    select(subject, locus, lineage, hit) %>%
    arrange(subject, lineage) %>%
    group_by(subject, locus, lineage) %>%
    summarise(hit = sum(hit)) %>%
    ungroup() %>%
    inner_join(rna, by = c("subject", "locus")) %>%
    filter(!is.na(rna)) %>%
    inner_join(distinct(dat_filtered, subject, locus), by = c("subject", "locus"))

estimate_df <- lm_df %>%
    split(.$locus) %>%
    map(~spread(., lineage, hit)) %>%
    map(~select(., -subject, -locus)) %>%
    lapply(function(x) coef(lm(rna ~ . - 1, data = x))) %>%
    map(. %>% as.data.frame %>% as_tibble) %>%
    map(~rownames_to_column(., "lineage")) %>%
    bind_rows(.id = "locus") %>%
    select(locus, lineage, y = 3)

fitted_df <- lm_df %>%
    split(.$locus) %>%
    map(~spread(., lineage, hit)) %>%
    map(~mutate(., f = paste("rna ~", paste(names(.)[-(1:3)], collapse = "+"), "- 1"),
		y = fitted(lm(as.formula(f))))) %>%
    map(~select(., subject, locus, rna, y)) %>%
    bind_rows() %>%
    left_join(genos, by = c("subject", "locus")) %>%
    group_by(subject, locus, rna, y) %>%
    summarise(lineage = paste(lineage, collapse = "/")) %>%
    ungroup() %>%
    select(subject, locus, lineage, rna, y)

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

