devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(scales)


loci <- gencode_hla$gene_name

hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% loci) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>% 
    left_join(geuvadis_info, by = c("subject" = "ena_id")) %>%
    select(subject = name, locus, tpm) %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))),
           source = "HLA-personalized")

hla_geuvadis <-  
    "../geuvadis_reanalysis/data/quantifications/peer/published/geuvadis_fpkms.csv" %>%
    read_csv() %>%
    filter(subject %in% hlapers$subject) %>%
    select(subject, gencode_hla_v12$gene_id) %>%
    gather(gene_id, fpkm, -subject) %>%
    left_join(gencode_hla_v12, by = "gene_id") %>%
    select(subject, locus = gene_name, fpkm) %>%
    mutate(locus = factor(locus, levels = loci)) %>%
    complete(subject, locus, fill = list(fpkm = 0)) %>%
    mutate(locus = reorder(locus, fpkm, median),
           locus = factor(locus, levels = rev(levels(locus))),
           source = "Original Geuvadis")

p1 <- ggplot(hlapers, aes(locus, tpm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "TPM")

p2 <- ggplot(hla_geuvadis, aes(locus, fpkm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "FPKM")

tiff("./plots/S1_fig.tiff", width = 4, height = 6, units = "in", res = 300)
plot_grid(p1, p2, ncol = 1)
dev.off()
