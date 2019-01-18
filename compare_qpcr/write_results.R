library(tidyverse)

pcr <- 
    bind_rows(read_tsv("HLA_A_formatted.tsv"),
              read_tsv("HLA_B_formatted.tsv"),
              read_tsv("HLA_C_formatted.tsv") %>% select(lineage, mrna)) %>%
    group_by(lineage) %>%
    filter(n() >= 10) %>%
    ungroup() %>%
    mutate(locus = sub("^([ABC])\\*\\d+$", "\\1", lineage)) %>%
    select(locus, lineage, mrna)

pcr_a <- filter(pcr, locus == "A")
pcr_b <- filter(pcr, locus == "B")
pcr_c <- filter(pcr, locus == "C")

pcr_pvals_a <- 
    pairwise.wilcox.test(pcr_a$mrna, pcr_a$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals_b <- 
    pairwise.wilcox.test(pcr_b$mrna, pcr_b$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals_c <- 
    pairwise.wilcox.test(pcr_c$mrna, pcr_c$lineage, p.adjust.method = "fdr") %>% 
    broom::tidy() %>%
    mutate(a1 = pmin(group1, group2), a2 = pmax(group1, group2)) %>%
    select(a1, a2, p = p.value)

pcr_pvals <- bind_rows(pcr_pvals_a, pcr_pvals_b, pcr_pvals_c)

pcr_medians <- pcr %>%
    group_by(lineage) %>%
    summarise(median = median(mrna)) %>%
    ungroup()

pcr_pairs <- tibble(locus = pcr$locus, a1 = pcr$lineage, a2 = pcr$lineage) %>%
    group_by(locus) %>%
    expand(a1, a2) %>%
    filter(a1 < a2) %>%
    ungroup() %>%
    left_join(pcr_medians, by = c("a1" = "lineage")) %>%
    left_join(pcr_medians, by = c("a2" = "lineage")) %>%
    mutate(sign = sign(median.x - median.y)) %>%
    select(a1, a2, sign)


hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(locus = sub("HLA-", "", locus),
           allele = sub("IMGT_", "", allele),
           lineage = sub("^([^:]+).+$", "\\1", allele)) %>%
    select(locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() >= 10) %>%
    ungroup()

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
    mutate(sign = sign(median.x - median.y)) %>%
    select(a1, a2, sign)

pvals_comp <- 
    left_join(pcr_pvals, hlapers_pvals, by = c("a1", "a2"), suffix = c("_qpcr", "_hlapers")) %>%
    mutate(signif_qpcr = ifelse(p_qpcr <= 0.05, 1, 0), 
           signif_hlapers = ifelse(p_hlapers <= 0.05, 1, 0)) %>%
    left_join(pcr_pairs, by = c("a1", "a2")) %>%
    left_join(hlapers_pairs, by = c("a1", "a2"), suffix = c("_qpcr", "_hlapers"))

write_tsv(pvals_comp, "results.tsv")

pvals_comp %>%
    filter(signif_qpcr == 1, signif_hlapers == 1) %>%
    mutate(same_dir = sign_qpcr == sign_hlapers) %>%
    arrange(desc(same_dir), a1, a2) %>%
    mutate(i = seq_len(n())) %>%
    select(i, everything()) %>%
    write_tsv("results_significant.tsv")

