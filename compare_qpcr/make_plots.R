library(tidyverse)

make_plot <- function(plot_df) {
    ggplot(plot_df, aes(reorder(to_order, rna, FUN = median), rna)) +
        ggbeeswarm::geom_quasirandom(aes(color = hom), size = .25, alpha = 1/3,
                                     method = "smiley", varwidth = TRUE, 
                                     show.legend = FALSE) +
        stat_summary(fun.y = median, geom = "point", shape = "\U2014", size = 4,
                     color = "cornflowerblue") +
        scale_x_discrete(labels = function(x) sub("_.+$", "", x)) +
        scale_y_continuous(labels = function(x) ifelse(x > 999, scales::comma(x), x)) +
        scale_color_manual(values = c("0" = "black", "1" = "green")) +
        facet_wrap(~locus+method, ncol = 1, scales = "free", strip.position = "left") +
        theme_bw() +
        theme(text = element_text(size = 8, family = "Times"), 
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              strip.text = element_text(face = "bold")) +
        labs(x = NULL, y = "mRNA")
}


pcr <- read_tsv("./data/hla_carrington_data_filteredN10.tsv") %>%
    group_by(subject, locus) %>%
    mutate(hom = as.integer(n() == 2 & n_distinct(lineage) == 1)) %>%
    ungroup() %>%
    arrange(subject, locus, lineage) %>%
    select(locus, lineage, hom, rna)

hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(locus = sub("HLA-", "", locus),
           allele = sub("IMGT_", "", allele),
           lineage = sub("^([^:]+).+$", "\\1", allele)) %>%
    group_by(lineage) %>%
    filter(n_distinct(subject) >= 10) %>%
    group_by(subject, locus) %>%
    mutate(hom = as.integer(n() == 2 & n_distinct(lineage) == 1)) %>%
    ungroup() %>%
    select(locus, lineage, hom, rna = tpm) 

plot_df <- list(qPCR = pcr, RNAseq = hlapers) %>%
    bind_rows(.id = "method") %>%
    unite(to_order, c("lineage", "method"), sep = "_", remove = FALSE) %>%
    mutate(hom = factor(hom)) %>%
    filter(lineage %in% unique(hlapers$lineage) & lineage %in% unique(pcr$lineage))


png("./FigSA.png", width = 8.5, height = 14, units = "cm", res = 300)
make_plot(plot_df)
dev.off()
