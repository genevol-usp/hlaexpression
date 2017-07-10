devtools::load_all("~/genomicRutils")
devtools::load_all("~/hlatools")
library(tidyverse)

# Functions
readLOG <- function(log) {
  x <- readLines(log)
  q <- x[grep("\\[quant\\] processed", x)]
  n <- as.integer(gsub(",", "", sub("^.+reads, ([0-9,]+).+$", "\\1", q)))
}

plot_lineages <- function(data) {
  ggplot(data, aes(x = reorder(lineage, expression, FUN = mean, na.rm = TRUE), 
                   y = expression)) +
    ggbeeswarm::geom_quasirandom(aes(color = factor(hom)), 
                                 method = "smiley", varwidth = TRUE,
                                 size = .75, show.legend = FALSE) +
    scale_color_manual(values = c("0" = "grey25", "1" = "green")) +
    stat_summary(fun.y = mean, geom = "point", shape = "\U2014", size = 9) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = NULL, y = "qPCR expression") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          legend.text = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 16))
}

# Data
samples <- readLines("./nci_samples.txt")

library_df <- read_tsv("./library_sizes.tsv") %>%
  rename(total = libsize)

log_files <- file.path("./quantifications_2/log", paste0(samples, ".log"))
names(log_files) <- samples

aligned <- map_df(log_files, readLOG) %>%
  gather(subject, aligned)

aligned_df <- inner_join(library_df, aligned, by = "subject") %>%
  gather(reads, value, -subject) %>%
  arrange(subject, reads)

nci_allele <- read_tsv("./nci_expression.tsv")

nci_gene <- distinct(nci_allele, subject, locus, mRNA, c_surface) %>%
  filter(subject %in% samples)

nci_lineage <- nci_allele %>%
  mutate(lineage = sub("^([^:]+).+$", "\\1", allele)) %>%
  select(subject, locus, lineage, expression = mRNA) %>%
  group_by(subject, locus) %>%
  mutate(hom = ifelse(length(unique(lineage)) == 1, 1L, 0L)) %>%
  filter(!is.na(expression)) %>%
  group_by(lineage) %>%
  filter(n() > 5) %>%
  ungroup() %>%
  mutate(fit = lm(expression ~ lineage - 1)$fitted.values)
      
rnaseq_allele <- read_tsv("./quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
  mutate(allele = hla_trimnames(gsub("IMGT_|_s\\d", "", allele), 3)) %>%
  select(subject, locus, allele, tpm)

rnaseq_gene <- rnaseq_allele %>%
  group_by(subject, locus) %>%
  summarize(tpm = sum(tpm))

rnaseq_lineage <- rnaseq_allele %>%
  mutate(lineage = hla_trimnames(allele, 1)) %>%
  group_by(subject, locus) %>%
  mutate(hom = ifelse(length(unique(allele)) == 1, 1L, 0L)) %>%
  select(subject, locus, lineage, hom, expression = tpm)

gene_df <- inner_join(rnaseq_gene, nci_gene, by = c("subject", "locus"))

corr <- gene_df %>%
  group_by(locus) %>%
  summarize(cor = round(cor(mRNA, tpm), digits = 3),
            rsq = round(summary(lm(tpm ~ mRNA))$adj.r.squared, digits = 3),
            x1 = min(mRNA) + 0.01,
            y1 = max(tpm) - 10,
            y2 = max(tpm) - 110) %>%
  ungroup() %>%
  mutate(label1 = paste("italic(r) ==", cor),
         label2 = paste("italic(r[adj]^{2}) == ", rsq))

hlac_df <- gene_df %>%
  filter(locus == "HLA-C", !is.na(c_surface)) %>%
  select(subject, rnaseq = tpm, qPCR = mRNA, c_surface) %>%
  gather(method, expression, rnaseq, qPCR) %>%
  select(subject, method, expression, c_surface) %>%
  arrange(subject, method)
  
corr_c <- hlac_df %>%
  group_by(method) %>%
  summarize(cor = round(cor(expression, c_surface), digits = 3),
            rsq = round(summary(lm(c_surface ~ expression))$adj.r.squared, digits = 3),
            x1 = min(expression) + 0.01,
            y1 = 500,
            y2 = 450) %>%
  ungroup() %>%
  mutate(label1 = paste("italic(r) ==", cor),
         label2 = paste("italic(r[adj]^{2}) == ", rsq))

rnaseq_hk <- read_tsv("./quantifications_2/housekeeping_norm_expression.tsv") %>%
  group_by(subject, locus) %>%
  summarize(tpm = sum(tpm), 
            est_counts_ACTBnorm = sum(est_counts_ACTBnorm),
            est_counts_B2Mnorm = sum(est_counts_B2Mnorm),
            est_counts_GAPDHnorm = sum(est_counts_GAPDHnorm))

cov_files <- file.path("./alignments", samples, "imgt.coverage")
names(cov_files) <- samples

covs_df <- 
  plyr::ldply(cov_files, 
              . %>% read_tsv(col_names = c("allele", "pos", "cov")) %>%
                filter(grepl("IMGT_(A|B|C)", allele)),
              .id = "subject", .parallel = TRUE) %>% as_tibble() %>%
  mutate(allele = sub("IMGT_", "", allele),
         locus = sub("^([^*]+).+$", "HLA-\\1", allele)) %>%
  select(subject, locus, allele, pos, cov)

haps <-
  covs_df %>%
  distinct(subject, locus, allele) %>%
  group_by(subject, locus) %>%
  mutate(h = paste0("A", 1:n())) %>%
  ungroup()

covs_df <-
  left_join(covs_df, haps, by = c("subject", "locus", "allele")) %>%
  select(subject, locus, h, allele, pos, cov)

covs_gene <- 
  covs_df %>%
  group_by(subject, locus, pos) %>%
  summarize(cov = sum(cov)) %>%
  ungroup()

covs_gene_100to300 <- 
  filter(covs_gene, pos >= 100 & pos <= 300) %>%
  group_by(subject, locus) %>%
  summarize(cov = mean(cov)) %>%
  ungroup()

qpcr <- select(nci_gene, subject, locus, qPCR = mRNA)

kallisto <- 
  read_tsv("./quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("HLA-A", "HLA-B", "HLA-C")) %>%
  group_by(subject, locus) %>%
  summarize(RNAseq_counts = sum(est_counts), RNAseq_TPM = sum(tpm)) %>%
  ungroup()

cov_expression_df <- 
  left_join(covs_gene_100to300, kallisto, by = c("subject", "locus")) %>%
  left_join(qpcr, by = c("subject", "locus")) %>%
  gather(method, value, RNAseq_counts, RNAseq_TPM, qPCR) %>%
  mutate(method = factor(method, levels = c("qPCR", "RNAseq_TPM", "RNAseq_counts")))


# plots
png("../plots/library_sizes.png", width = 12, height = 5, units = "in", res = 300)
ggplot(aligned_df) +
  geom_point(aes(x = reorder(subject, value, FUN = "max"), y = value,
                color = reads), size = 1) +
  ggsci::scale_color_aaas() +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 250e6, 25e6)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  labs(x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
dev.off()

png("../plots/seq_vs_pcr.png", width=12, height=4, units="in", res=300)
ggplot(gene_df, aes(mRNA, tpm)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~locus, scales = "free") +
  geom_text(data = corr, aes(x = x1, y = y1, label = label1), 
            parse = TRUE, hjust = 0, size = 6) +
  geom_text(data = corr, aes(x = x1, y = y2, label = label2), 
            parse = TRUE, hjust = 0, size = 6) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 16)) +
  labs(x = "qPCR", y = "RNAseq (TPM)")
dev.off()

png("../plots/ab_vs_rna.png", width=10, height=4, units="in", res=300)
ggplot(hlac_df, aes(expression, c_surface)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~method, scales = "free") +
  geom_text(data = corr_c, aes(x = x1, y = y1, label = label1), 
            parse = TRUE, hjust = 0, size = 6) +
  geom_text(data = corr_c, aes(x = x1, y = y2, label = label2), 
            parse = TRUE, hjust = 0, size = 6) +
  theme_bw() +
  labs(x = "RNA", y = "Surface Protein") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 16))
dev.off()

png("../plots/rnaseq_lineages.png", width=12, height=6, units="in", res=300)
plot_lineages(rnaseq_lineage)
dev.off()

png("../plots/nci_lineages.png", width=12, height=6, units="in", res=300)
plot_lineages(nci_lineage)
dev.off()

png("../plots/tpm_vs_B2MnormCounts.png", width = 12, height = 6, units = "in", res = 300)
ggplot(rnaseq_hk, aes(tpm, est_counts_B2Mnorm)) +
  geom_point(size = 2, alpha = .5) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~locus, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18)) +
  labs(x = "TPM", y = "B2M-normalized counts") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = 6,
                        label.x.npc = .3, label.y.npc = .7) 
dev.off()

png("../plots/tpm_vs_ACTBnormCounts.png", width = 12, height = 6, units = "in", res = 300)
ggplot(rnaseq_hk, aes(tpm, est_counts_ACTBnorm)) +
  geom_point(size = 2, alpha = .5) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~locus, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18)) +
  labs(x = "TPM", y = "ACTB-normalized counts") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = 6,
                        label.x.npc = .7, label.y.npc = .05) 
dev.off()

png("../plots/tpm_vs_GAPDHnormCounts.png", width = 12, height = 6, units = "in", res = 300)
ggplot(rnaseq_hk, aes(tpm, est_counts_GAPDHnorm)) +
  geom_point(size = 2, alpha = .5) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~locus, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18)) +
  labs(x = "TPM", y = "GAPDH-normalized counts") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = 6,
                        label.x.npc = .7, label.y.npc = .05) 
dev.off()

png("../plots/covs_measure.png", width = 10, height = 5, units = "in", res = 300)
ggplot(data = filter(covs_df, subject == '66K00634')) +
  geom_line(data = filter(covs_gene, subject == '66K00634'), 
            aes(pos, cov, group = locus)) +
  geom_line(aes(pos, cov, group = allele, color = h), show.legend = FALSE) +
  ggrepel::geom_label_repel(data = filter(covs_df, subject == '66K00634', pos == 600L),
                            aes(pos, cov, group = allele, color = h, label = allele),
                            show.legend = FALSE) +
  ggsci::scale_color_aaas() +
  facet_wrap(~locus, ncol = 3, scales = "free") +
  labs(x = "position", y = "coverage") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))
dev.off()

png("../plots/covs_expression.png", width = 8, height = 8, units = "in", res = 300)
ggplot(cov_expression_df, aes(value, cov)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
  facet_wrap(~locus+method, scales = "free", ncol = 3) +
  labs(x = "expression", y = "average coverage (positions 100bp to 300bp)") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = 4,
                        label.x.npc = .3, label.y.npc = .7) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12))
dev.off()
