devtools::load_all("~/hlaseqlib/")
library(tidyverse)
library(cowplot)

loci <- readLines("../../../imgt_index_v2/imgt_loci.txt") %>%
    paste0("HLA-", .)

hla_personalized <- 
    "../../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
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

gencode12_hla <-
    "~/gencode_data/gencode.v12.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% loci)

phen <- read_csv("./peer/published/geuvadis_fpkms.csv")

hla_geuvadis <- phen %>%
    filter(subject %in% hla_personalized$subject) %>%
    select(subject, gencode12_hla$gene_id) %>%
    gather(gene_id, fpkm, -subject) %>%
    left_join(gencode12_hla, by = "gene_id") %>%
    select(subject, locus = gene_name, fpkm) %>%
    mutate(locus = factor(locus, levels = loci)) %>%
    complete(subject, locus, fill = list(fpkm = 0)) %>%
    mutate(locus = reorder(locus, fpkm, median),
           locus = factor(locus, levels = rev(levels(locus))),
           source = "Original Geuvadis")

p1 <- ggplot(hla_personalized, aes(locus, tpm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")

p2 <- ggplot(hla_geuvadis, aes(locus, fpkm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12)) +
    labs(y = "FPKM")

png("./expression_boxplot.png", width = 6, height = 6, units = "in", res = 300)
plot_grid(p1, p2, ncol = 1)
dev.off()
    
