devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)

# functions
calc_ase <- function(counts) min(counts)/sum(counts)

plot_lower <- function(data, mapping, ...) {
    
    ggplot(data = data, mapping = mapping) +
	geom_point(size = .5) +
	geom_smooth(method = lm) 
}

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

make_data <- function(locus1, locus2, df) {
  
    cis <- select(df, subject, locus1, locus2) %>%
        drop_na() %>%
        rename(gene1 = !!locus1, gene2 = !!locus2)

    trans <- cis %>% 
	    group_by(subject) %>%
	    mutate_at(vars(gene2), rev) %>%
	    ungroup()

    gene <- cis %>%
        group_by(subject) %>%
        summarize_at(vars(gene1, gene2), sum) %>%
        ungroup()

    bind_rows(list("Gene-level" = gene, 
                   "Within haplotypes" = cis, 
                   "Between haplotypes" = trans), .id = "level") %>%
	    mutate(pair = paste(locus1, "vs", locus2),
	        level = factor(level, levels = c("Gene-level", 
	                                            "Within haplotypes", 
	                                            "Between haplotypes"))) %>%
        select(pair, everything())
        
}

# data
hla_genes <- sort(gencode_hla$gene_name)

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

## TPM
kallisto_imgt_tpm <- 
    read_tsv("./kallisto/imgt/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup()

kallisto_pri_tpm <- 
    read_tsv("./kallisto/pri/quantifications/processed_imgt_quants.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm) %>%
    arrange(subject, locus)

star_imgt <- 
    read_tsv("./star/imgt/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele))

star_imgt_tpm <-
    star_imgt %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup()

star_imgt_with_nonclassical <- 
    read_tsv("./star/imgt/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(grepl("HLA", locus)) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

star_pri_tpm <- 
    read_tsv("./star/pri/quantifications/processed_imgt_quants.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm)

## PCA
kallisto_imgt_pca <-
    "../qtls/kallisto/pca/imgt/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

kallisto_pri_pca <- 
    "../qtls/kallisto/pca/pri/phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_imgt_pca <-
    "../qtls/star/imgt/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_pri_pca <- 
    "../qtls/star/pri/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

quant_data <- 
    left_join(kallisto_imgt_tpm, star_imgt_tpm, by = c("subject", "locus"),
              suffix = c(".kallisto.imgt", ".star.imgt")) %>%
    left_join(kallisto_pri_tpm, by = c("subject", "locus")) %>%
    rename(tpm.kallisto.pri = tpm) %>%
    left_join(star_pri_tpm, by = c("subject", "locus")) %>%
    rename(tpm.star.pri = tpm) %>%
    left_join(kallisto_imgt_pca, by = c("subject", "locus")) %>%
    rename(pca.kallisto.imgt = resid) %>%
    left_join(star_imgt_pca, by = c("subject", "locus")) %>%
    rename(pca.star.imgt = resid) %>%
    left_join(kallisto_pri_pca, by = c("subject", "locus")) %>%
    rename(pca.kallisto.pri = resid) %>%
    left_join(star_pri_pca, by = c("subject", "locus")) %>%
    rename(pca.star.pri = resid)

star_tpm_df <- quant_data %>%
    select(subject, locus, tpm.star.imgt, tpm.star.pri) %>%
    gather(index, tpm, 3:4) %>%
    mutate(index = sub("^tpm\\.star\\.", "", index)) %>%
    arrange(subject, locus, index)

pag3f <- mutate(pag, allele = hla_trimnames(allele, 3))

ase_df <- star_imgt %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarize(ase = calc_ase(est_counts)) %>%
    ungroup()
  
ase_error <- star_imgt %>%
    select(subject, locus, allele) %>%
    mutate(locus = sub("^HLA-", "", locus),
           allele = hla_trimnames(allele)) %>%
    calc_genotyping_accuracy(pag3f, by_locus = FALSE) %>%
    group_by(subject, locus) %>%
    summarize(error = sum(!correct)) %>%
    ungroup() %>%
    mutate(locus = paste0("HLA-", locus)) %>%
    inner_join(ase_df, by = c("subject", "locus"))

hla_and_transAct_genes <- gencode_chr_gene %>%
    filter(gene_name %in% c(gencode_hla$gene_name, "CIITA"))

classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "../expression/star/imgt/quantifications_expressed50%.bed" %>%
    read_tsv() %>%
    filter(gid %in% classII_and_CIITA$gene_id) %>%
    select(gid, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, tpm, -1) %>%
    left_join(hla_and_transAct_genes, by = c("gid" = "gene_id")) %>%
    select(subject, locus = gene_name, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    select(subject, locus, tpm, CIITA) %>%
    arrange(subject, locus)

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

pca_star_df <-
    sprintf("../qtls/star/imgt/1-phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
    setNames(pcs) %>%
    map_df(. %>% read_tsv(progress = FALSE) %>% select(gid, matches("^HG|^NA")), 
	   .id = "PC") %>%
    gather(subject, value, -(PC:gid))

pca_star_hla <-
    inner_join(pca_star_df, hla_and_transAct_genes, by = c("gid" = "gene_id")) %>%
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

concordant_haps_I <- 
    read_tsv("./star/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
    select(subject, hap)

concordant_haps_II <- 
    read_tsv("./star/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
    select(subject, hap)

tpm_by_allele_I <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_I, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_allele_II <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(!locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_II, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_allele <- 
    full_join(tpm_by_allele_I, tpm_by_allele_II, by = c("subject", "hap"))

# plots
png("./plots/star_vs_kallisto_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.star.imgt", "tpm.kallisto.imgt") +
    labs(x = "TPM (STAR-Salmon)", y = "TPM (kallisto)")
dev.off()

png("./plots/star_vs_kallisto_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "pca.star.imgt", "pca.kallisto.imgt") +
    labs(x = "PCA-corrected estimate (STAR-Salmon)", y = "PCA-corrected estimate (kallisto)")
dev.off()

png("./plots/star_imgt_vs_pri_TPM.png", height = 6, width = 10, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.star.imgt", "tpm.star.pri") +
    labs(x = "TPM (HLA personalized index)", y = "TPM (Ref transcriptome)")
dev.off()

png("./plots/star_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "pca.star.imgt", "pca.star.pri") +
    labs(x = "PCA-corrected estimate (HLA personalized index)", 
         y = "PCA-corrected estimate (Ref transcriptome)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.kallisto.imgt", "tpm.kallisto.pri") +
    labs(x = "TPM (HLA personalized index)", y = "TPM (Ref transcriptome)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "pca.kallisto.imgt", "pca.kallisto.pri") +
    labs(x = "PCA-corrected estimate (HLA personalized index)", 
         y = "PCA-corrected estimate (Ref transcriptome)")
dev.off()

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 200)
ggplot(star_tpm_df, aes(tpm, fill = index)) +
    geom_density(alpha = 1/2) +
    scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
    ggthemes::scale_fill_colorblind() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, margin = margin(t = 10))) +
    facet_wrap(~locus, scales = "free")
dev.off()

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

png("./plots/hlacorrelations.png", width = 6, height = 6, units = "in", res = 200)
pairs_hla_k <-
    star_imgt_tpm %>%
    select(-est_counts) %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>%
    select(-subject) %>%
    ggpairs(lower = list(continuous = plot_lower), upper = list()) + 
    theme_bw() +
    theme(title = element_text(size = 14),
          axis.text.x = element_text(angle = 90))

print(pairs_hla_k, left = .3, bottom = .3)
dev.off()

png("./plots/trans_activ_corrs.png", width = 10, height = 3.5, units = "in", res = 200)
cor_df <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$tpm, .$CIITA),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(r = round(r, digits = 3))

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

png("./plots/within_vs_between_haps.png", height = 8, width = 8, units = "in", res = 200)

plotphase <- function(phase_df) {
    
    phase_cor_df <- phase_df %>%
        group_by(pair, level) %>%
        summarize(r = cor(gene1, gene2),
                  x = min(gene1),
                  y = max(gene2)) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 2))
    
    ggplot(phase_df, aes(gene1, gene2)) +
        geom_point(size = 1) +
        geom_smooth(method = lm, se = FALSE) + 
        geom_text(data = phase_cor_df,
                  aes(x, y, label = paste("r =", r)),
                  hjust = "inward", vjust = "inward", size = 4) +
        facet_grid(pair~level, scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_blank(),
              strip.text = element_text(size = 10))
}

phase_data <- 
    tribble(
        ~locus1, ~locus2,
        "A"    , "B",
        "A"    , "C",
        "B"    , "C",
        "DQA1" , "DQB1",
        "DQA1" , "DRB1") %>%
    pmap_df(make_data, tpm_by_allele)

phase_list <- filter(phase_data, level != "Gene-level") %>%
    split(.$pair)

p1 <- plotphase(phase_list[[1]])
p2 <- plotphase(phase_list[[2]]) + theme(strip.text.x = element_blank())    
p3 <- plotphase(phase_list[[3]]) + theme(strip.text.x = element_blank())   
p4 <- plotphase(phase_list[[4]]) + theme(strip.text.x = element_blank())   
p5 <- plotphase(phase_list[[5]]) + theme(strip.text.x = element_blank())

plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, .9, .9, .9, .9))

dev.off()


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

png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 200)
ggplot(star_imgt_with_nonclassical, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()
