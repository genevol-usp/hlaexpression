devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggrepel)


read_conditional_mhc <- function(file, mhc_genes) {
    
    qtltools_res <- 
        read_qtltools(file) %>%
        inner_join(mhc_genes, by = c("phen_id" = "gene_id")) %>%
        mutate(tss = ifelse(strand == "+", phen_from, phen_to),
               dist_to_tss = var_from - tss,
               pval = -log10(bwd_pval)) %>%
        select(gene = gene_name, tss, strand, rank, var_id, 
               var_from, dist, dist_to_tss, pval, best = bwd_best, 
                   signif = bwd_signif) 
    
    fix_rank <- qtltools_res %>%
        filter(best == 1) %>%
        arrange(gene, desc(pval), rank) %>%
        group_by(gene) %>%
        mutate(new_rank = 1:n() - 1L) %>%
        ungroup() %>%
        select(gene, rank, new_rank)
    
    out_df <- left_join(qtltools_res, fix_rank, by = c("gene", "rank")) %>%
        select(gene, tss, strand, rank = new_rank, var_id, var_from, dist,
               dist_to_tss, pval, best, signif) %>%
        group_by(gene, var_id) %>%
        filter(pval == max(pval)) %>%
        ungroup() %>%
        mutate(rank = as.character(rank))
    
    out_df
}


mhc_coords <- gencode_hla %>% 
    summarise(start = min(start) - 1e6L, end = max(end) + 1e6L)

mhc_genes <- gencode_pri_gene %>%
    filter(start >= mhc_coords$start, end <= mhc_coords$end) %>%
    select(gene_id, gene_name)

hlapers_qtls <-
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/2-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

ref_qtls <-
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/reference/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

qtls_df <- 
    list("HLApers" = hlapers_qtls, "Ref Transcriptome" = ref_qtls) %>%
    bind_rows(.id = "index") 

hla_qtls <- hlapers_qtls %>%
    mutate(gene = factor(gene, levels = unique(gene)),
           colorize = ifelse(gene %in% gencode_hla$gene_name, 1L, 0L))


plot_qtls_1 <- ggplot() +
    geom_point(data = filter(hla_qtls, colorize == 0L), 
               aes(var_from, pval), 
               color = "grey", alpha = .1, size = .5) +
    geom_point(data = filter(hla_qtls, colorize == 1L),
               aes(var_from, pval, color = gene),
               alpha = .25, size = .5) +
    coord_cartesian(xlim = c(29.2e6L, 33.8e6L)) +
    ggsci::scale_color_npg(labels = function(x) sub("HLA-", "", x)) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    guides(color = guide_legend(keyheight = .5, override.aes = list(alpha = 1, size = 2.5))) +
    theme_bw() +
    theme(text = element_text(family = "Times", size = 9)) +
    labs(x = "position (Mb)",
         y = expression(paste("-log"[10], italic(Pvalue)))) 


plot_qtls_df <- qtls_df %>%
    filter(gene %in% gencode_hla$gene_name) %>%
    mutate(gene = sub("HLA-", "", gene),
           gene = factor(gene, levels = sub("HLA-", "", gencode_hla$gene_name)))

ori_df <- plot_qtls_df %>%
    distinct(index, gene, strand) %>%
    mutate(xend = ifelse(strand == "+", 1e5L, -1e5L)) %>%
    select(index, gene, xend)

best_df <- filter(plot_qtls_df, best == 1L)

plot_qtls_2 <-  ggplot(data = plot_qtls_df, aes(dist_to_tss, pval)) +
    geom_blank(data = plot_qtls_df %>% 
                   group_by(gene, index) %>%
                   slice(which.max(pval)) %>% 
                   mutate(max_pval = pval + 1L),
               aes(y = max_pval)) +
    coord_cartesian(xlim = c(-7e5, +7e5)) +
    geom_vline(aes(xintercept = 0), color = "grey", size = 1.5) + 
    geom_point(data = filter(plot_qtls_df, signif == 0), 
               aes(dist_to_tss, pval),
               size = .5, color = "grey", show.legend = FALSE) +
    geom_segment(data = ori_df,
                 aes(x = 0, xend = xend, y = -2.5, yend = -2.5),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed", ends = "last"),
                 size = 1, color = "darkgreen", alpha = .5) +
    geom_point(data = filter(plot_qtls_df, signif == 1L), 
               aes(dist_to_tss, pval, color = rank), 
               size = .5) +
    geom_point(data = best_df, 
               aes(dist_to_tss, pval, color = rank), size = 2) +
    geom_point(data = best_df, 
               aes(dist_to_tss, pval), 
               shape = 1, size = 2, color = "black", stroke = 1) +
    geom_label_repel(data = best_df,
                    aes(dist_to_tss, pval, label = var_id),
                    label.size = NA, label.padding = 0.05, alpha = 0.4, seed = 1,
                    size = 2.5, family = "Times", fontface = "bold",
                    direction = "both",
                    show.legend = FALSE) +
    geom_label_repel(data = best_df,
                     aes(dist_to_tss, pval, label = var_id),
                     label.size = NA, label.padding = 0.05, fill = NA, seed = 1,
                     size = 2.5, family = "Times", fontface = "bold",
                     direction = "both",
                     show.legend = FALSE) +
    scale_color_manual(values = c("0" = "#8491B4B2",
                                  "1" = "#DC0000B2",
                                  "2" = "gold3",
                                  "3" = "#7E6148FF",
                                  "4" = "#009E73")) +
    scale_x_continuous(breaks = seq(-5e5L, 5e5L, by = 2.5e5L),
                       labels = function(x) x/1e6L) +
    facet_grid(gene~index, scales = "free_y") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    theme_bw() +
    theme(text = element_text(family = "Times", size = 9)) +
    labs(x = "distance from TSS (Mb)", 
         y = expression(paste("-log"[10], italic(Pvalue))))


tiff("./plots/Fig5.tiff", width = 12, height = 18, units = "cm", res = 300)
plot_grid(NULL, plot_qtls_1, plot_qtls_2, ncol = 1,
          rel_heights = c(.025, .225, 1),
          labels = c("A", "", "B"), label_size = 9, label_fontfamily = "Times",
          vjust = 1.4, hjust = -0.7)
dev.off()


#### non-classical

gencode_hla_nclassical <- gencode_chr_gene %>% 
    filter(gene_name %in% c("HLA-E", "HLA-F", "HLA-G", "HLA-H"))

plot_qtls_df <- qtls_df %>%
    filter(gene %in% gencode_hla_nclassical$gene_name) %>%
    mutate(gene = sub("HLA-", "", gene),
           gene = factor(gene, levels = sub("HLA-", "", gencode_hla_nclassical$gene_name)))

ori_df <- plot_qtls_df %>%
    distinct(index, gene, strand) %>%
    mutate(xend = ifelse(strand == "+", 1e5L, -1e5L)) %>%
    select(index, gene, xend)

best_df <- filter(plot_qtls_df, best == 1L)

png("./plots/eqtls_nonclassic.png", width = 12, height = 9, units = "cm", res = 150)
ggplot(data = plot_qtls_df, aes(dist_to_tss, pval)) +
    geom_blank(data = plot_qtls_df %>% 
                   group_by(gene, index) %>%
                   slice(which.max(pval)) %>% 
                   mutate(max_pval = pval + 1L),
               aes(y = max_pval)) +
    coord_cartesian(xlim = c(-1e6, +1e6)) +
    geom_vline(aes(xintercept = 0), color = "grey", size = 1.5) + 
    geom_point(data = filter(plot_qtls_df, signif == 0), 
               aes(dist_to_tss, pval),
               size = .5, color = "grey", show.legend = FALSE) +
    geom_segment(data = ori_df,
                 aes(x = 0, xend = xend, y = -2.5, yend = -2.5),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed", ends = "last"),
                 size = 1, color = "darkgreen", alpha = .5) +
    geom_point(data = filter(plot_qtls_df, signif == 1L), 
               aes(dist_to_tss, pval, color = rank), 
               size = .5) +
    geom_point(data = best_df, 
               aes(dist_to_tss, pval, color = rank), size = 2) +
    geom_point(data = best_df, 
               aes(dist_to_tss, pval), 
               shape = 1, size = 2, color = "black", stroke = 1) +
    geom_label_repel(data = best_df,
                     aes(dist_to_tss, pval, label = var_id),
                     label.size = NA, label.padding = 0.05, alpha = 0.4, seed = 1,
                     size = 2.5, family = "Times", fontface = "bold",
                     direction = "both",
                     show.legend = FALSE) +
    geom_label_repel(data = best_df,
                     aes(dist_to_tss, pval, label = var_id),
                     label.size = NA, label.padding = 0.05, fill = NA, seed = 1,
                     size = 2.5, family = "Times", fontface = "bold",
                     direction = "both",
                     show.legend = FALSE) +
    scale_color_manual(values = c("0" = "#8491B4B2",
                                  "1" = "#DC0000B2",
                                  "2" = "gold3",
                                  "3" = "#7E6148FF",
                                  "4" = "#009E73")) +
    scale_x_continuous(breaks = seq(-1e6L, 1e6L, by = 2.5e5L),
                       labels = function(x) x/1e6L) +
    scale_y_continuous(breaks = seq(0, 30, by = 5)) +
    facet_grid(gene~index, scales = "free_y") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    theme_bw() +
    theme(text = element_text(family = "Times", size = 9)) +
    labs(x = "distance from TSS (Mb)", 
         y = expression(paste("-log"[10], italic(Pvalue))))
dev.off()
