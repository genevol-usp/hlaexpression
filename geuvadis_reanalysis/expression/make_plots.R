devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
library(scales)

# data
hla_genes <- gencode_hla$gene_name

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

# Boxplot
star_all_imgt_loci <- 
    read_tsv("./star/supplemented/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(grepl("HLA", locus)) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 200)
ggplot(star_all_imgt_loci, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()


# Scatterplots of pipeline comparisons

scatter_plot_cors <- function(df, x_var, y_var) {
    
    cor_df <- df %>%
        group_by(locus) %>%
        do(data.frame(r = cor(.[[x_var]], .[[y_var]]),
                      x = min(.[[x_var]]),
                      y = max(.[[y_var]]))) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 3))
        
    ggplot(df, aes_string(x_var, y_var)) +
        geom_abline() +
        geom_point(size = .8) +
        geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
                  hjust = "inward", vjust = "inward", size = 5) +
        facet_wrap(~locus, scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16))
}


## STAR vs kallisto
kallisto_imgt <- 
    read_tsv("./kallisto/supplemented/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

star_imgt <- 
    read_tsv("./star/supplemented/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

star_kallisto_df <-
    left_join(star_imgt, kallisto_imgt, by = c("subject", "locus"), 
	      suffix = c(".star", ".kallisto"))

png("./plots/star_vs_kallisto_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(star_kallisto_df, "tpm.star", "tpm.kallisto") +
    labs(x = "TPM (STAR-Salmon)", y = "TPM (kallisto)")
dev.off()


## STAR vs kallisto (PCA)
kallisto_imgt_pca <-
    "../qtls/kallisto/supplemented/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_imgt_pca <-
    "../qtls/star/supplemented/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_kallisto_pca_df <-
    left_join(star_imgt_pca, kallisto_imgt_pca, by = c("subject", "locus"), 
	      suffix = c(".star", ".kallisto"))

png("./plots/star_vs_kallisto_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(star_kallisto_pca_df, "resid.star", "resid.kallisto") +
    labs(x = "PCA-corrected estimate (STAR-Salmon)", y = "PCA-corrected estimate (kallisto)")
dev.off()


## STAR: supplemented vs reference transcriptome
star_pri <- 
    read_tsv("./star/transcriptome/quantifications/processed_imgt_quants.tsv") %>%
    select(subject, locus, tpm)

star_imgt_vs_pri_df <- 
    left_join(star_imgt, star_pri, by = c("subject", "locus"), 
	      suffix = c(".imgt", ".ref"))

png("./plots/star_imgt_vs_pri_TPM.png", height = 6, width = 10, units = "in", res = 200)
scatter_plot_cors(star_imgt_vs_pri_df, "tpm.imgt", "tpm.ref") +
    labs(x = "TPM (HLA personalized index)", y = "TPM (Ref transcriptome)")
dev.off()


## STAR: supplemented vs reference transcriptome (PCA)
star_pri_pca <- 
    "../qtls/star/transcriptome/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_imgt_vs_pri_pca_df <- 
    left_join(star_imgt_pca, star_pri_pca, by = c("subject", "locus"),
	      suffix = c(".imgt", ".ref"))

png("./plots/star_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(star_imgt_vs_pri_pca_df, "resid.imgt", "resid.ref") +
    labs(x = "PCA-corrected estimate (HLA personalized index)", 
         y = "PCA-corrected estimate (Ref transcriptome)")
dev.off()


# TPM distributions

tpm_distribution_df <- star_imgt_vs_pri_df %>%
    gather(index, tpm, 3:4) %>%
    mutate(index = sub("^tpm\\.", "", index)) %>%
    arrange(subject, locus, index)

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 200)
ggplot(tpm_distribution_df, aes(tpm, fill = index)) +
    geom_density(alpha = 1/2) +
    scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
    scale_fill_manual(values = c(imgt = "#8491B4B2", ref = "#DC0000B2"),
		       labels = c(imgt = "HLA-supplemented", ref = "Ref transcriptome")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, margin = margin(t = 10))) +
    facet_wrap(~locus, scales = "free")
dev.off()


# ASE
calc_ase <- function(counts) min(counts)/sum(counts)

pag3f <- mutate(pag, allele = hla_trimnames(allele, 3))

df_for_ase <- 
    "./star/supplemented/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    ungroup() %>%
    left_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts) %>%
    mutate(allele = gsub("IMGT_", "", allele),
	   allele = hla_trimnames(allele))

ase_df <- df_for_ase %>% 
    group_by(subject, locus) %>%
    summarize(ase = calc_ase(est_counts)) %>%
    ungroup() 
  
ase_error <- df_for_ase %>%
    select(-est_counts) %>%
    mutate(locus = sub("^HLA-", "", locus)) %>%
    calc_genotyping_accuracy(pag3f, by_locus = FALSE) %>%
    group_by(subject, locus) %>%
    summarize(error = sum(!correct)) %>%
    ungroup() %>%
    mutate(locus = paste0("HLA-", locus)) %>%
    inner_join(ase_df, by = c("subject", "locus"))

png("./plots/ase.png", width = 8, height = 5, units = "in", res = 200)
ggplot(ase_df, aes(locus, ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1, alpha = 1/2) +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))
dev.off()

png("./plots/ase_genot_errors.png", width = 8, height = 4, units = "in", res = 200)
ggplot(ase_error, aes(locus, ase, color = factor(error))) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1.5, alpha = 1/2) +
    coord_cartesian(ylim = c(0, 0.5)) +
    scale_color_manual(values = c("0" = "grey", "1" = "blue", "2" = "red")) +
    labs(x = "", color = "genotyping errors") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))
dev.off()

png("./plots/ase_histogram.png", width = 8, height = 4, units = "in", res = 200)
ggplot(ase_df, aes(ase)) +
  geom_density(fill = "grey35", color = NA) +
  facet_wrap(~locus) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
dev.off()


# Correlations

## CRD
gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, start, end, strand)

crd <- "../../LCL_ALL.chr6.subset.txt.gz" %>%
    read_delim(col_names = FALSE, delim = " ") %>%
    select(-X2, -X4, -X6, -X8, -X9, -X11) %>%
    filter(between(X3, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 1e6),
           between(X7, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 5e5)) %>%
    mutate(index = (X1 + X5)/2L,
           pos = (X3 + X7)/2L,
           y = X5 - X1,
           r2 = X10^2) %>%
    group_by(index) %>%
    mutate(pos_label = min(pos)) %>%
    ungroup()

gene_pos <- gencode_hla_v19 %>%
    mutate(gene_name = sub("HLA-", "", gene_name),
           pos = ifelse(strand == "+", start, end),
           closest = crd$index[map_dbl(pos, ~which.min(abs(. - crd$pos_label)))])

crd_plot <- ggplot() +
    geom_point(data = crd, 
               aes(x = index, y = y, alpha = r2), 
               color = "blue", size = .5) +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_x_continuous(breaks = crd$index[seq(1, nrow(crd), 4e4)],
                       labels = round(crd$pos_label[seq(1, nrow(crd), 4e4)]/1e6, 1)) +
    geom_hline(yintercept = -15, size = 2, color = "grey", alpha = 1/2) +
    geom_segment(data = filter(gene_pos, strand == "+"),
                 aes(x = closest, xend = closest+7.5, y = -15, yend = -15),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last"),
                 alpha = 1/2) +
    geom_segment(data = filter(gene_pos, strand == "-"),
                 aes(x = closest, xend = closest-7.5, y = -15, yend = -15),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last"),
                 color = "red", alpha = 1/2) +
    geom_text(data = filter(gene_pos, strand == "+"),
              aes(x = closest, y = -5, label = gene_name),
              size = 3, hjust = 0) +
    ggrepel::geom_text_repel(data = filter(gene_pos, strand == "-"),
                             aes(x = closest, y = -25, label = gene_name),
                             size = 3, color = "red") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = c(.95, .8)) +
    labs(x = "Genomic position (Mb)", alpha = expression(~r^2)) +
    guides(alpha = guide_legend(override.aes = list(size = 2.5)))
 

## Global correlations
global_cors <- star_imgt %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>%
    select(-subject) %>% ggcorr(label = TRUE)

## Within vs between haplotypes

make_data <- function(locus1, locus2, df) {
    
    cis <- select(df, subject, locus1, locus2) %>%
        drop_na() %>%
        rename(gene1 = !!locus1, gene2 = !!locus2)
    
    trans <- cis %>% 
        group_by(subject) %>%
        mutate_at(vars(gene2), rev) %>%
        ungroup()
    
    bind_rows(list("Within haplotypes" = cis, 
                   "Between haplotypes" = trans), .id = "level") %>%
        mutate(pair = paste(locus1, "vs", locus2),
               level = factor(level, levels = c("Within haplotypes", "Between haplotypes"))) %>%
        select(pair, everything())
}

plotphase <- function(phase_df) {
    
    phase_cor_df <- phase_df %>%
        group_by(pair, level) %>%
        summarize(r = cor(gene1, gene2),
                  x = min(gene1),
                  y = max(gene2)) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 2))
    
    ggplot(phase_df, aes(gene1, gene2)) +
        geom_point(size = .7) +
        geom_smooth(method = lm, se = FALSE) + 
        scale_x_continuous(breaks = scales::pretty_breaks(2)) +
        scale_y_continuous(breaks = scales::pretty_breaks(2)) +
        geom_text(data = phase_cor_df,
                  aes(x, y, label = paste("r =", r)),
                  hjust = "inward", vjust = "inward", size = 4) +
        facet_grid(pair~level, scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_blank(),
              strip.text = element_text(size = 8))
}


haps_data <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_rsid.tsv") %>%
    filter(!grepl("/", hla_allele)) %>%
    mutate(locus = sub("HLA-", "", locus),
           tpm = as.numeric(tpm)) %>%
    select(subject, hap, locus, tpm) %>% 
    spread(locus, tpm)

phase_data <- 
    tribble(
        ~locus1, ~locus2,
        "A"    , "B",
        "A"    , "C",
        "B"    , "C",
        "DQA1" , "DQB1",
        "DQA1" , "DRB1") %>%
    pmap_df(make_data, haps_data)

phase_list <- filter(phase_data, level != "Gene-level") %>%
    split(.$pair)


p1 <- plotphase(phase_list[[1]])
p2 <- plotphase(phase_list[[2]]) + theme(strip.text.x = element_blank())    
p3 <- plotphase(phase_list[[3]]) + theme(strip.text.x = element_blank())   
p4 <- plotphase(phase_list[[4]]) + theme(strip.text.x = element_blank())   
p5 <- plotphase(phase_list[[5]]) + theme(strip.text.x = element_blank())

within_vs_between <- 
    plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, .9, .9, .9, .9))

plot_AB <- plot_grid(crd_plot, global_cors, ncol = 1, rel_heights = c(.55, .45))

png("./plots/correlations.png", width = 10, height = 6.5, units = "in", res = 300)
plot_grid(plot_AB, within_vs_between, ncol = 2, rel_widths = c(.55, .45))
dev.off()

# Transactivator

classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "../expression/star/supplemented/quantifications_expressed50%.bed" %>%
    read_tsv() %>%
    inner_join(classII_and_CIITA, by = c("gid" = "gene_id")) %>%
    select(gid, locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, tpm, -gid, -locus) %>%
    select(subject, locus, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    select(subject, locus, tpm, CIITA) %>%
    arrange(subject, locus)

cor_df <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$tpm, .$CIITA),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(r = round(r, digits = 3))

png("./plots/trans_activ_corrs.png", width = 10, height = 3, units = "in", res = 200)
ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = 1) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16)) + 
    facet_wrap(~locus, nrow = 1, scales = "free") +
    labs(x = NULL)
dev.off()


# Correlation decrease as number of PCs increase

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

pca_star_df <-
    sprintf("../qtls/star/supplemented/1-phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
    setNames(pcs) %>%
    map_df(. %>% read_tsv(progress = FALSE) %>% select(gid, matches("^HG|^NA")), 
	   .id = "PC") %>%
    gather(subject, value, -(PC:gid))

hla_and_ciita <- gencode_chr_gene %>%
    filter(gene_name %in% c(hla_genes, "CIITA"))

pca_star_hla <-
    inner_join(pca_star_df, hla_and_ciita, by = c("gid" = "gene_id")) %>%
    select(PC, subject, gene_name, value) %>%
    mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
    spread(gene_name, value) 

cors_pca_star <- pca_star_hla %>%
    group_by(PC) %>%
    summarize(AxB = cor(A, B), 
              AxC = cor(A, C), 
              BxC = cor(B, C),
              DQA1xDQB1 = cor(DQA1, DQB1), 
              DQA1xDRB1 = cor(DQA1, DRB1), 
              DQB1xDRB1 = cor(DQB1, DRB1),  
              DQA1xCIITA = cor(DQA1, CIITA),
              DQB1xCIITA = cor(DQB1, CIITA), 
              DRB1xCIITA = cor(DRB1, CIITA)) %>%
    gather(gene_pair, correlation, -1) %>%
    mutate(PC = as.integer(PC)) %>%
    arrange(PC, gene_pair)

png("./plots/correlation_decrease.png", width = 10, height = 5, units = "in", res = 200)
ggplot(cors_pca_star, aes(PC, correlation)) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = pcs) +
    facet_wrap(~gene_pair) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.position = "top",	
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          axis.text.x = element_text(angle = 90)) +
    labs(x = "Number of PCs")
dev.off()
