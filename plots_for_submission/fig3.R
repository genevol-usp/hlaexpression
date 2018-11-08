devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)

hla_genes <- gencode_hla$gene_name

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

#imgt_loci <- readLines("../imgt_index/imgt_loci.txt") %>%
#    paste0("HLA-", .)

mapping_imgt <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    #filter(locus %in% imgt_loci) %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    group_by(locus) %>%
    filter(mean(tpm) >= 100) %>%
    ungroup() %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

#tiff("./plots/Fig3.tiff", width = 5, height = 3, units = "in", res = 300)
tiff("./plots/Fig3.tiff", width = 4, height = 3, units = "in", res = 300)
ggplot(mapping_imgt, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = comma) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Arial"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "TPM")
dev.off()