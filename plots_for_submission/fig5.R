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
        mutate(gene = reorder(gene, tss),
               rank = as.character(rank))
    
    out_df
}

plot_qtls_index <- function(qtl_df) {
    
   
    
   
}

hla_qtl_genes <- gencode_hla %>% 
    filter(!gene_name %in% c("HLA-DRA", "HLA-DPA1")) %>%
    pull(gene_name)

mhc_coords <- gencode_hla %>% 
    filter(gene_name %in% hla_qtl_genes) %>%
    summarise(start = min(start) - 1e6L, end = max(end) + 1e6L)

mhc_genes <- gencode_pri_gene %>%
    filter(start >= mhc_coords$start, end <= mhc_coords$end) %>%
    select(gene_id, gene_name)

hla_qtls <-
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/2-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

ref_qtls <-
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/reference/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

qtls_df <- 
    list("HLA-personalized" = hla_qtls, "Ref Transcriptome" = ref_qtls) %>%
    bind_rows(.id = "index") %>%
    mutate(gene = reorder(gene, tss))

hla_qtls <- hla_qtls %>%
    mutate(gene = reorder(gene, tss),
           colorize = ifelse(gene %in% hla_qtl_genes, 1L, 0L))

plot_qtls_1 <- ggplot() +
    geom_point(data = filter(hla_qtls, colorize == 0L), 
               aes(var_from, pval), 
               color = "grey", alpha = .1, size = .5) +
    geom_point(data = filter(hla_qtls, colorize == 1L),
               aes(var_from, pval, color = gene),
               alpha = .25, size = .5) +
    coord_cartesian(xlim = c(29.2e6L, 33.8e6L)) +
    ggsci::scale_color_npg() +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme(text = element_text(size = 9, family = "Arial"),
          axis.text = element_text(size = 8, family = "Arial"),
          axis.title.y = element_text(hjust = 0)) +
    guides(color = guide_legend(keyheight = .5, override.aes = list(alpha = 1, size = 3))) +
    labs(x = "position (Mb)",
         y = expression(paste("-log"[10], italic(Pvalue)))) 


plot_qtls_df <- qtls_df %>%
    filter(gene %in% hla_qtl_genes) 

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
    geom_text_repel(data = best_df,
                    aes(dist_to_tss, pval, label = var_id),
                    fontface = "bold", size = 3, 
                    segment.color = "gray35",
                    nudge_x = 1.75e5, force = 1,
                    show.legend = FALSE) +
    scale_color_manual(values = c("0" = "#8491B4B2",
                                  "1" = "#DC0000B2",
                                  "2" = "gold3",
                                  "3" = "#7E6148FF",
                                  "4" = "#009E73")) +
    scale_x_continuous(breaks = seq(-5e5L, 5e5L, by = 2.5e5L),
                       labels = function(x) x/1e6L) +
    theme_bw() +
    theme(text = element_text(size = 10, family = "Arial"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey")) +
    facet_grid(gene~index, scales = "free_y") +
    labs(x = "distance from TSS (Mb)", 
         y = expression(paste("-log"[10], italic(Pvalue)))) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))



tiff("./plots/Fig5.tiff",  width = 6.25, height = 7.25, units = "in", res = 300)
plot_grid(NULL, plot_qtls_1, plot_qtls_2, ncol = 1,
          rel_heights = c(.025, .225, 1),
          labels = c("", "A", "B"), label_size = 12, vjust = 1, hjust = 0)
dev.off()
