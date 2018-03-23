devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(ggplot2)

# Functions
plot_lineages <- function(data) {
    ggplot(data, aes(x = reorder(lineage, expression, FUN = mean, na.rm = TRUE), 
		     y = expression)) +
	ggbeeswarm::geom_quasirandom(aes(color = factor(hom)), 
				     method = "smiley", varwidth = TRUE,
				     size = .75, show.legend = FALSE) +
	scale_color_manual(values = c("0" = "grey25", "1" = "green")) +
	stat_summary(fun.y = mean, geom = "point", shape = "\U2014", size = 9) +
	facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
	labs(x = NULL, y = "mRNA") +
	theme_bw() +
	theme(axis.title = element_text(size = 16),
	      axis.text.y = element_text(size = 14),
	      axis.text.x = element_text(size = 12),
	      legend.text = element_text(size = 14),
	      strip.text = element_text(face = "bold", size = 16))
}

# Data
samples <- readLines("./data/nci_sample_ids.txt")

nci_allele <- read_tsv("./data/nci_expression.tsv")

nci_gene <- distinct(nci_allele, subject, locus, mRNA, c_surface) %>%
    filter(subject %in% samples)

nci_lineage <- nci_allele %>%
    mutate(lineage = hla_trimnames(allele, 1)) %>%
    select(subject, locus, lineage, expression = mRNA) %>%
    group_by(subject, locus) %>%
    mutate(hom = ifelse(length(unique(lineage)) == 1, 1L, 0L)) %>%
    filter(!is.na(expression)) %>%
    group_by(lineage) %>%
    filter(n() > 5) %>%
    ungroup() %>%
    mutate(est_allele_exp = lm(expression ~ lineage - 1)$fitted.values)
      
rnaseq_allele <- 
    read_tsv("./expression/star/quantifications_2/processed_quant.tsv") %>%
    filter(grepl("HLA", locus)) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    select(subject, locus, allele, tpm)

rnaseq_gene <- rnaseq_allele %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

rnaseq_lineage <- rnaseq_allele %>%
    mutate(lineage = hla_trimnames(allele, 1)) %>%
    group_by(subject, locus) %>%
    mutate(hom = ifelse(length(unique(allele)) == 1, 1L, 0L)) %>%
    select(subject, locus, lineage, hom, expression = tpm) %>%
    ungroup()

common <-
    inner_join(select(rnaseq_lineage, subject, locus, lineage),
               select(nci_lineage, subject, locus, lineage),
               by = c("subject", "locus", "lineage")) %>%
    distinct()

gene_df <- inner_join(rnaseq_gene, nci_gene, by = c("subject", "locus"))

hlac_df <- gene_df %>%
  filter(locus == "HLA-C", !is.na(c_surface)) %>%
  select(subject, rnaseq = tpm, qPCR = mRNA, c_surface) %>%
  gather(method, expression, rnaseq, qPCR) %>%
  select(subject, method, expression, c_surface) %>%
  arrange(subject, method)
  
qpcr <- select(nci_gene, subject, locus, qPCR = mRNA)

# plots
png("./plots/seq_vs_pcr.png", width=12, height=4, units="in", res=300)
ggplot(gene_df, aes(mRNA, tpm)) +
  geom_point(size = 2) +
  facet_wrap(~locus, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 16)) +
  labs(x = "qPCR", y = "RNAseq (TPM)") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
			formula = y ~ x, parse = TRUE, size = 6)
dev.off()

png("./plots/ab_vs_rna.png", width=10, height=4, units="in", res=300)
ggplot(hlac_df, aes(expression, c_surface)) +
  geom_point(size = 2) +
  facet_wrap(~method, scales = "free") +
  theme_bw() +
  labs(x = "RNA", y = "Surface Protein") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 16)) +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
			  formula = y ~ x, parse = TRUE, size = 6)
dev.off()

png("./plots/rnaseq_lineages.png", width=12, height=6, units="in", res=300)
rnaseq_lineage %>%
    inner_join(common) %>%
    plot_lineages() + labs(y = "TPM")
dev.off()

png("./plots/nci_lineages.png", width=12, height=6, units="in", res=300)
nci_lineage %>%
    inner_join(common) %>%
    plot_lineages() + labs(y = "mRNA (qPCR)")
dev.off()

png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 200)
ggplot(rnaseq_gene, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 90, 
				     hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()
