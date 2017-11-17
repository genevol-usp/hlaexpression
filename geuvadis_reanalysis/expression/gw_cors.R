devtools::load_all("~/hlaseqlib")
library(tidyverse)

bed <- read_tsv("../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur.bed.gz")

tpm_df <- bed %>%
    select(gid, matches("^HG|^NA")) %>%
    gather(subject, tpm, -gid)

qt4_genes <- tpm_df %>%
    group_by(gid) %>%
    summarize(tpm = mean(tpm)) %>%
    ungroup() %>%
    filter(tpm >= 100) %>%
    pull(gid)

tibble(gene1 = qt4_genes, gene2 = qt4_genes) %>%
    expand(gene1, gene2) %>%
    filter(gene1 != gene2) %>%
    group_by(gene1) %>%
    sample_n(1) %>%
    ungroup() %>%
    left_join(tpm_df, by = c("gene1" = "gid")) %>%
    left_join(tpm_df, by = c("subject", "gene2" = "gid"), suffix = c("1", "2")) %>%
    group_by(gene1, gene2) %>%
    summarize(r = cor(tpm1, tpm2)) %>%
    ungroup() %>%
    summarize(mean_r = mean(r))

