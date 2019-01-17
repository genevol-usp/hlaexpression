library(tidyverse)

reorder_within <- function(x, by, within, fun = median, sep = "___", ...) {
    new_x <- paste(x, within, sep = sep)
    reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
    reg <- paste0(sep, ".+$")
    ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

make_plot <- function(plot_df) {
    ggplot(plot_df, aes(reorder_within(lineage, mrna, method, FUN = median), mrna)) +
        ggbeeswarm::geom_quasirandom(color = ggsci::pal_npg()(6)[6], size = .5,
                                     method = "smiley", varwidth = TRUE,
                                     show.legend = FALSE) +
        geom_boxplot(outlier.shape = NA, fill = NA, alpha = .5) +
        scale_x_reordered() +
        facet_wrap(~method, ncol = 1, scales = "free") +
        theme_bw() +
        theme(text = element_text(size = 8, family = "Times"), 
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              strip.text = element_text(face = "bold")) +
        labs(x = NULL, y = "mRNA")
}


pcr <- 
    bind_rows(read_tsv("HLA_A_formatted.tsv"),
              read_tsv("HLA_B_formatted.tsv"),
              read_tsv("HLA_C_formatted.tsv") %>% select(lineage, mrna)) %>%
    mutate(locus = sub("^([ABC])\\*\\d+$", "\\1", lineage)) %>%
    select(locus, lineage, mrna)

hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
    mutate(locus = sub("HLA-", "", locus),
           allele = sub("IMGT_", "", allele),
           lineage = sub("^([^:]+).+$", "\\1", allele)) %>%
    select(locus, lineage, mrna = tpm) 

plot_df <- list(qPCR = pcr, RNAseq = hlapers) %>%
    bind_rows(.id = "method")

plot_df_a <- filter(plot_df, locus == "A")
plot_df_b <- filter(plot_df, locus == "B")
plot_df_c <- filter(plot_df, locus == "C")


tiff("./FigSA.png", width = 10, height = 10, units = "cm", res = 300)
make_plot(plot_df_a)
dev.off()

tiff("./FigSB.png", width = 10, height = 10, units = "cm", res = 300)
make_plot(plot_df_b)
dev.off()

tiff("./FigSC.png", width = 10, height = 10, units = "cm", res = 300)
make_plot(plot_df_c)
dev.off()
