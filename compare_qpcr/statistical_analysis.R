library(tidyverse)

pcr <- read_tsv("./data/hla_carrington_data_filteredN10.tsv")

pcr_a <- filter(pcr, locus == "A")
pcr_b <- filter(pcr, locus == "B")
pcr_c <- filter(pcr, locus == "C")

pcr_pvals_a <- 
    pairwise.wilcox.test(pcr_a$rna, pcr_a$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals_b <- 
    pairwise.wilcox.test(pcr_b$rna, pcr_b$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals_c <- 
    pairwise.wilcox.test(pcr_c$rna, pcr_c$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals <- bind_rows(pcr_pvals_a, pcr_pvals_b, pcr_pvals_c)

pcr_medians <- pcr %>%
    group_by(lineage) %>%
    summarise(median = median(rna)) %>%
    ungroup()

pcr_pairs <- tibble(locus = pcr$locus, a1 = pcr$lineage, a2 = pcr$lineage) %>%
    group_by(locus) %>%
    expand(a1, a2) %>%
    filter(a1 < a2) %>%
    ungroup() %>%
    left_join(pcr_medians, by = c("a1" = "lineage")) %>%
    left_join(pcr_medians, by = c("a2" = "lineage")) %>%
    mutate(direction = sign(median.x - median.y)) %>%
    select(a1, a2, direction)


hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(locus = sub("HLA-", "", locus),
           allele = sub("IMGT_", "", allele),
           lineage = sub("^([^:]+).+$", "\\1", allele)) %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    mutate(nLineage = n_distinct(subject)) %>%
    group_by(subject, locus) %>%
    filter(all(nLineage >= 10)) %>%
    ungroup() %>%
    select(-nLineage) %>%
    arrange(subject, locus, lineage)

hlapers_a <- filter(hlapers, locus == "A")
hlapers_b <- filter(hlapers, locus == "B")
hlapers_c <- filter(hlapers, locus == "C")

hlapers_pvals_a <- 
    pairwise.wilcox.test(hlapers_a$tpm, hlapers_a$lineage, p.adjust.method = "fdr") %>%
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

hlapers_pvals_b <- 
    pairwise.wilcox.test(hlapers_b$tpm, hlapers_b$lineage, p.adjust.method = "fdr") %>%
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

hlapers_pvals_c <- 
    pairwise.wilcox.test(hlapers_c$tpm, hlapers_c$lineage, p.adjust.method = "fdr") %>%
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

hlapers_pvals <- bind_rows(hlapers_pvals_a, hlapers_pvals_b, hlapers_pvals_c)

hlapers_medians <- hlapers %>%
    group_by(lineage) %>%
    summarise(median = median(tpm)) %>%
    ungroup()

hlapers_pairs <- 
    tibble(locus = hlapers$locus, a1 = hlapers$lineage, a2 = hlapers$lineage) %>%
    group_by(locus) %>%
    expand(a1, a2) %>%
    filter(a1 < a2) %>%
    ungroup() %>%
    left_join(hlapers_medians, by = c("a1" = "lineage")) %>%
    left_join(hlapers_medians, by = c("a2" = "lineage")) %>%
    mutate(direction = sign(median.x - median.y)) %>%
    select(a1, a2, direction)

pvals_comp <- 
    left_join(pcr_pvals, hlapers_pvals, 
	      by = c("a1", "a2"), suffix = c("_qpcr", "_hlapers")) %>%
    mutate(signif_qpcr = ifelse(p_qpcr <= 0.05, 1, 0), 
           signif_hlapers = ifelse(p_hlapers <= 0.05, 1, 0)) %>%
    left_join(pcr_pairs, by = c("a1", "a2")) %>%
    left_join(hlapers_pairs, by = c("a1", "a2"), suffix = c("_qpcr", "_hlapers"))

significant <- pvals_comp %>%
    filter(signif_qpcr == 1, signif_hlapers == 1) %>%
    select(-starts_with("signif"))

write_tsv(significant, "results_significant.tsv")


summary_df <- significant %>%
    mutate(locus = sub("^([^\\*])\\*.+$", "\\1", a1)) %>%
    group_by(locus) %>%
    summarise(same_dir = sum(direction_qpcr == direction_hlapers),
	      n_significant = n()) %>%
    ungroup()

write_tsv(summary_df, "summary_results.tsv")
