devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)


dist_ref <- read_tsv("../imgt_index/distance_to_ref/distances_to_reference.tsv")

hlapers_dist <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(dist_ref, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarise(dist = mean(dist), hlapers = sum(tpm)) %>%
    ungroup()

ref_transc <-
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    select(subject, locus, ref = tpm)

express_df <- left_join(hlapers_dist, ref_transc, by = c("subject", "locus")) %>%
    mutate(locus = sub("HLA-", "", locus),
           locus = factor(locus, levels = sub("HLA-", "", gencode_hla$gene_name)))

cor_df <- express_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$hlapers, .$ref),
                  p = cor(.$hlapers, .$ref, method = "spearman"),
                  x = min(.$hlapers),
                  y = max(.$ref))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 2)) %>%
    mutate(label = paste("r == ", r))


tiff("./plots/Fig4.tiff", width = 12, height = 8, units = "cm", res = 300)
ggplot(express_df, aes(hlapers, ref)) +
    geom_abline() +
    geom_point(aes(color = dist), size = .6) +
    scale_x_continuous(labels = comma, breaks = pretty_breaks(2)) +
    scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
    scale_color_gradient(low = "white", high = "darkred") +
    geom_text(data = cor_df, aes(x, y, label = label), parse = TRUE, 
              hjust = "inward", vjust = "inward", size = 2.5, 
              color = "white") +
    facet_wrap(~locus, scales = "free") +
    guides(color = guide_colorbar(barheight = 4)) +
    theme_dark() +
    theme(text = element_text(size = 9, family = "Times"), 
          strip.text = element_text(face = "bold")) +
    labs(x = "HLApers", y = "Ref Transcriptome", color = "Divergence\nfrom Ref")
dev.off()