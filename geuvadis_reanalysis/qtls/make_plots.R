devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggpmisc)
library(ggplot2)

# genotype PCA 
make_pca_plot <- function(PC_x, PC_y) {
  ggplot(pcs, aes_string(PC_x, PC_y)) +
    geom_point(aes(color = pop)) +
    scale_color_manual(values = c(GBR = "darkolivegreen4", FIN = "dodgerblue2",
                                  CEU = "red", TSI = "purple", YRI = "goldenrod")) +
    theme_bw()
}

geuvadis_pops <- select(geuvadis_info, subject = name, pop)

pcs <- 
  read_delim("./pca_genotypes/eur.pca", delim = " ") %>%
  mutate(SampleID = sub("^.+_(PC\\d+)$", "\\1", SampleID)) %>%
  filter(SampleID %in% paste0("PC", 1:100)) %>%
  gather(subject, value, -1) %>%
  spread(SampleID, value) %>%
  mutate_at(vars(-subject), function(x) x/sqrt(sum(x^2))) %>%
  inner_join(geuvadis_pops, by = "subject")

png("./plots/genotype_pca.png", width = 10, height = 5, units = "in", res = 300)
p1 <- make_pca_plot("PC1", "PC2")
p2 <- make_pca_plot("PC2", "PC3")
p3 <- make_pca_plot("PC3", "PC4")
p4 <- make_pca_plot("PC4", "PC5")

leg <- get_legend(p1)

plot_grid(p1 + guides(color=FALSE), 
          p2 + guides(color=FALSE), 
          leg,
          p3 + guides(color=FALSE), 
          p4 + guides(color=FALSE), 
          nrow = 2, ncol = 3, rel_widths = c(3, 3, 1))
dev.off()

# Number of eQTLs per method of phenotype correction
egenes_pca <-
  sprintf("./qtls_kallisto/qtltools_correction/permutations/results/permutations_%d.significant.txt", seq(0, 100, 5)) %>%
  setNames(seq(0, 100, 5)) %>%
  parallel::mclapply(function(x) read_delim(x, delim = " ", col_names = FALSE), 
                     mc.cores = 21) %>%
  bind_rows(.id = "f") %>%
  count(f)

peer_perm_files <- 
  list.files("./qtls_kallisto/peer_correction/permutations/results", 
	     pattern = "permutations_\\d+\\.significant", full.names = TRUE)

names(peer_perm_files) <- 
  sub("^.+permutations_(\\d+)\\.significant\\.txt$", "\\1", peer_perm_files)

egenes_peer <-
  peer_perm_files %>%
  parallel::mclapply(function(x) read_delim(x, delim = " ", col_names = FALSE),
                     mc.cores = length(peer_perm_files)) %>%
  bind_rows(.id = "f") %>%
  count(f)

egenes_df <- list(PCA = egenes_pca, PEER = egenes_peer) %>%
  bind_rows(.id = "method") %>%
  mutate(f = as.integer(f)) %>%
  arrange(f, method)

png("./plots/pca_vs_peer.png", width = 8, height = 5, units = "in", res = 300)
ggplot(egenes_df, aes(f, n, color = method, group = method)) + 
  geom_point(size = 2.5) + 
  geom_line() +
  scale_x_continuous(breaks = sort(unique(egenes_df$f))) +
  ggsci::scale_color_npg() +
  theme_bw() +
  labs(x = "Number of PCs/factors", y = "Number of eGenes")
dev.off()

# Lineage and effects plot

## lineages
concordant <- 
  bind_rows(
    read_tsv("../expression/kallisto/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
      gather(locus, allele, A:C), 
    read_tsv("../expression/kallisto/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
      gather(locus, allele, DQA1:DRB1)) %>%
  arrange(subject, locus, hap) %>%
  mutate(concordant_phase = TRUE)

residuals_by_allele_best <- 
  read_tsv("../expression/kallisto/phase_hla_alleles/data/hla_allele_expression_bestpc.bed") %>%
  gather(subject, resid, -locus, -hap)

expression_data <-
  read_tsv("../expression/kallisto/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
  filter(rank == 0) %>%
  left_join(concordant, by = c("subject", "locus", "hap")) %>%
  left_join(residuals_by_allele_best, by = c("subject", "locus", "hap")) %>%
  select(subject, hap, locus, hla_allele, rank, resid,
                variant, variant_allele, concordant_phase)

exp_levels <-
  expression_data %>% 
  filter(!is.na(concordant_phase)) %>%
  group_by(locus, rank, variant, variant_allele) %>% 
  summarize(eQTL = mean(resid)) %>%
  group_by(locus, variant) %>%
  mutate(eQTL = ifelse(eQTL == max(eQTL), "High", "Low")) %>%
  ungroup()

k_levels <-
  expression_data %>%
  left_join(exp_levels, by = c("locus", "rank", "variant", "variant_allele")) %>%
  mutate(eQTL = ifelse(is.na(concordant_phase), NA_character_, eQTL)) %>%
  select(subject, locus, rank, hla_allele, eQTL, resid) %>%
  mutate(hla_allele = hla_trimnames(hla_allele, 1),
         eQTL = factor(eQTL, levels = c("High", "Low"))) %>%
  arrange(subject, locus, hla_allele)

## effects
gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_id, gene_name)

phen_best <-
  "./qtls_kallisto/qtltools_correction/phenotypes/phenotypes_eur_75.bed.gz" %>%
  read_tsv() %>%
  inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
  select(gene_name, HG00096:NA20828) %>%
  gather(subject, resid, -gene_name) %>%
  select(subject, locus = gene_name, resid)

best_eqtl_locus <-
  read_tsv("../expression/kallisto/phase_hla_alleles/best_eqtls.tsv") %>%
  filter(rank == 0L) %>%
  left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
  select(locus = gene_name, variant = var_id)

eqtl_info <- 
  read_tsv("../expression/kallisto/phase_hla_alleles/best_eqtl_snps.vcf", comment = "##") %>%
  select(-`#CHROM`, -POS, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  gather(subject, genotype, -(ID:ALT)) %>%
  inner_join(best_eqtl_locus, by = c("ID" = "variant")) %>%
  separate(genotype, c("h1", "h2"), convert = TRUE) %>%
  gather(hap, allele, h1:h2) %>%
  mutate(allele = ifelse(allele == 0, REF, ALT)) %>%
  arrange(subject, locus, ID, allele) %>%
  group_by(subject, locus, ID) %>%
  summarize(genotype = paste(sort(allele), collapse = "/")) %>%
  select(subject, locus, variant = ID, genotype) %>%
  ungroup() %>%
  unite(id, genotype, locus, sep = "_", remove = FALSE)

eqtls_expression_df <- left_join(eqtl_info, phen_best, by = c("subject", "locus"))

## plot
png("./plots/lineage_and_effects.png", width = 12, height = 12, units = "in", res = 300)
p1 <- 
  ggplot(k_levels, aes(x = reorder(hla_allele, resid, FUN = median, na.rm = TRUE), 
                       y = resid)) +
  geom_jitter(aes(color = eQTL), alpha = .5) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
  scale_color_manual(values = c("Low" = "royalblue", "High" = "red"), 
                     na.value = "grey") +
  scale_x_discrete(labels = function(x) gsub("\\*", "*\n", x)) +
  facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
  labs(x = " ", y = "TPM residualized by 60 PCs") +
  theme_bw() +
  theme(axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        strip.text = element_text(face = "bold", size = rel(1.2)),
        legend.position = "top") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

p2 <-
  ggplot(eqtls_expression_df, aes(reorder(id, resid, "mean"), resid)) +
  geom_jitter(width = .25, alpha = 1/2, size = .75) +
  geom_smooth(aes(group = 1), method = lm, se = FALSE) +
  scale_x_discrete(labels = function(x) paste0(sub("^([^_]+).+$", "\\1", x), "\n")) +
  scale_y_continuous(position = "right") +
  geom_text(data = distinct(eqtls_expression_df, locus, variant), 
            aes(x = 1.5, y = 3, label = variant)) +
  coord_cartesian(ylim = c(-3, 3.2)) +
  facet_wrap(~locus, ncol = 1, scales = "free") +
  labs(x = " ", y = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1)),
        strip.text = element_blank())

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

grid1 <- plot_grid(legend, NULL, p1, p2, ncol = 2, 
                   rel_widths = c(4, 1), rel_heights = c(.07, 1))

ggdraw(grid1) + 
  draw_label("HLA lineage", 0.44, 0.02, size = 16) +
  draw_label("1000G genotype", 0.88, 0.02, size = 16)
dev.off()

# eQTL landscape around TSS
conditional_pca <-
  read_qtltools("./qtls_kallisto/qtltools_correction/conditional_analysis/conditional_75_all.txt.gz") %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("phen_id" = "gene_id")) %>%
  mutate(dist_tss = ifelse(strand == "+", 
                           var_from - phen_from,
                           phen_to - var_from),
         nom_pval = -log10(bwd_pval)) %>%
  select(phen_id = gene_name, rank, var_id, dist_tss, nom_pval, 
         slope = bwd_slope, best = bwd_best, signif = bwd_signif) %>%
  group_by(phen_id, var_id, best) %>%
  filter(rank == min(rank)) %>%
  ungroup() %>%
  mutate(rank = factor(rank))

png("./plots/qtls_landscape.png", height = 12, width = 10, units = "in", res = 300)
ggplot() +
  geom_vline(xintercept = 0, color = "grey25", size = 2) + 
  geom_point(data = filter(conditional_pca, signif == 0), 
             aes(dist_tss, nom_pval),
             color = "grey", alpha = .1, show.legend = FALSE) +
  geom_point(data = filter(conditional_pca, signif == 1L), 
             aes(dist_tss, nom_pval, color = rank), 
             alpha = .5) +
  geom_point(data = filter(conditional_pca, best == 1L), 
             aes(dist_tss, nom_pval, color = rank), 
             size = 2) +
  geom_point(data = filter(conditional_pca, best == 1L), 
             aes(dist_tss, nom_pval), 
             shape = 1, size = 2, color = "black", stroke = 1.5) +
  geom_hline(data = conditional_pca %>% 
               group_by(phen_id) %>% 
               filter(signif == 0) %>% 
               summarise(thres = max(nom_pval)),
             aes(yintercept = thres), color = "black") +
  coord_cartesian(xlim = c(-1e6, +1e6)) +
  ggsci::scale_color_aaas() +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal() +
  facet_wrap(~phen_id, scales = "free_y", ncol = 1) +
  labs(x = "distance from TSS", 
       y = expression(paste("-log"[10], italic(Pvalue)))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
dev.off()
