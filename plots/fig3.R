devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)


hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(locus = sub("HLA-", "", locus),
           locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

tiff("./plots/Fig3.tiff", width = 10, height = 6, units = "cm", res = 300)
ggplot(hlapers, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = comma) +
    theme_bw() + 
    theme(text = element_text(family = "Times", size = 11)) +
    labs(x = NULL, y = "TPM")
dev.off()