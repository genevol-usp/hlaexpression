library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib")

dat <-
    "./data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz" %>%
    read_tsv(skip = 2)

hla <- filter(dat, Description %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(-1) %>%
    gather(tissue, tpm, -1) %>%
    rename(gene = Description) %>%
    group_by(tissue) %>%
    filter(all(tpm > 400)) %>%
    ungroup()

hlac <- hla %>%
    filter(gene == "HLA-C") %>%
    select(tissue, tpm_c = tpm)

hla_df <- left_join(hla, hlac, by = "tissue") %>%
    mutate(tpm_rel_c = tpm/tpm_c) %>%
    filter(gene != "HLA-C") %>%
    group_by(tissue) %>%
    mutate(ix = mean(tpm_rel_c)) %>%
    ungroup() %>%
    mutate(tissue = reorder(tissue, ix, mean)) %>%
    select(gene, tissue, tpm_rel_c)

png("./relative_tpm_gtex.png", width = 8, height = 5, units = "in", res = 300)
ggplot(hla_df, aes(tissue, tpm_rel_c, color = gene, group = gene)) +
    geom_line(size = 2) +
    ggsci::scale_color_npg() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    labs(x = NULL, y = "TPM relative to HLA-C") 
dev.off()

hla2 <- filter(dat, Description %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    select(-1) %>%
    gather(tissue, tpm, -1) %>%
    rename(gene = Description) 

hla_levels <- hla2 %>%
    group_by(tissue) %>%
    summarise(md = median(tpm)) %>%
    ungroup() %>%
    arrange(md)

hla2 <- left_join(hla2, hla_levels)

cols <- ggsci::pal_npg()(4) %>%
    .[c(1, 2, 4)] %>%
    setNames(c("HLA-A", "HLA-B", "HLA-C"))

png("./medians_gtex.png", width = 8, height = 5, units = "in", res = 300)
ggplot(hla2, aes(reorder(tissue, md, min), tpm, color = gene)) +
    geom_point(size = 3, alpha = .7) +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    labs(x = NULL, y = "median TPM")
dev.off()