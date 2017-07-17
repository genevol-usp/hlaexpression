devtools::load_all("/home/vitor/genomicRutils")
devtools::load_all("/home/vitor/hlatools")
library(tidyverse)
library(cowplot)
library(GGally)
library(ggpmisc)
library(ggplot2)

# library sizes
readLOG <- function(log) {
  x <- readLines(log)
  q <- x[grep("\\[quant\\] processed", x)]
  n <- as.integer(gsub(",", "", sub("^.+reads, ([0-9,]+).+$", "\\1", q)))
}

geuvadis_ids <- geuvadis_info %>%
  filter(kgp_phase3 == 1) %>%
  select(subject = ena_id, name)

library_size <- 
  read_delim("./library_size.txt", delim = " ") %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, total = n) %>%
  arrange(subject)

log_files <- 
  file.path("./quantifications_2/log", paste0(geuvadis_ids$subject, ".quant.log")) %>%
  setNames(geuvadis_ids$subject)

aligned <- lapply(log_files, readLOG)

aligned_df <- 
  tibble(subject = names(aligned), aligned = unlist(aligned)) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, aligned) %>%
  arrange(subject) %>%
  left_join(library_size, by = "subject") %>%
  gather(reads, n, aligned, total) %>%
  arrange(subject, reads)

png("./plots/library_sizes.png", height = 4, width = 10, units = "in", res = 150)
ggplot(aligned_df) +
  geom_line(aes(x = reorder(subject, n, FUN = "max"), y = n, 
                color = reads, group = reads), size = 1.1) +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()

# Comparison of expression estimates across methods

## kallisto vs hlaTX
hlatx <- 
 "~/hlaexpression/geuvadis_reanalysis/expression/hlaTX/results/results_TX.txt" %>%
  hla_readtx() %>%
  mutate(locus = sub("^([^_]+).*$", "\\1", allele)) %>%
  inner_join(geuvadis_ids) %>%
  select(subject = name, locus, est_counts) %>%
  group_by(subject, locus) %>%
  summarize(counts = sum(est_counts)) %>%
  ungroup()

expression_df <- 
  read_tsv("quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1")) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, allele, est_counts, tpm) %>%
  mutate(allele = hla_trimnames(gsub("IMGT_|_s\\d", "", allele), 3)) %>%
  arrange(subject, locus, allele)

kallisto_imgt <- expression_df %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts), tpm = sum(tpm))

kallisto_hlatx <- 
  inner_join(kallisto_imgt, hlatx, by = c("subject", "locus")) %>%
  select(subject, locus, kallisto = est_counts, hlatx = counts) %>%
  inner_join(library_size, by = "subject") %>%
  mutate(kallisto = kallisto/total, hlatx = hlatx/total) %>%
  select(-total)

png("./plots/kallisto_vs_hlatx.png", height = 6, width = 10, units = "in", res = 150)
ggplot(kallisto_hlatx, aes(kallisto, hlatx)) +
  geom_abline(color = "grey") +
  geom_point(alpha = 1/2) +
  scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
  scale_y_continuous(breaks = function(y) scales::pretty_breaks(3)(y)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1.3,
                                   margin = margin(t = 7)),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) +
  labs(x = "kallisto", y = "GEM-based method") +
  facet_wrap(~locus, scales = "free") +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
               formula = y ~ x, parse = TRUE, size = 5) 
dev.off()

## kallisto vs Geuvadis
gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")))

gencode12 <- 
  get_gencode_coords("~/gencode_data/gencode.v12.annotation.gtf.gz", feature = "gene") %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(locus = gene_name, Gene_Symbol = gene_id)

geuvadis <- 
  "~/hlaexpression/geuvadis_reanalysis/data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz" %>%
  read_tsv() %>%
  inner_join(gencode12, by = "Gene_Symbol") %>%
  select(locus, HG00096:NA20828) %>%
  gather(subject, fpkm, HG00096:NA20828) %>%
  filter(subject %in% geuvadis_ids$name) %>%
  select(subject, locus, fpkm) %>%
  arrange(subject, locus)

kallisto_fpkms <- 
  read_tsv("./quantifications_2/gene_fpkms.tsv") %>%
  left_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, fpkm)

fpkm_df <- 
  inner_join(kallisto_fpkms, geuvadis, by = c("subject", "locus")) %>%
  select(subject, locus, kallisto = fpkm.x, geuvadis = fpkm.y)

png("./plots/kallisto_vs_geuvadis.png", width = 10, height = 6, units = "in", res = 300)
ggplot(fpkm_df, aes(kallisto, geuvadis)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~locus, scales = "free") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                        formula = y ~ x, parse = TRUE, size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) +
  labs(x = "FPKM (kallisto)", y = "FPKM (GEUVADIS)")
dev.off()

## kallisto different indices
kallisto_chr <- 
  read_tsv("./quantifications_CHR/processed_quant.tsv") %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, tpm)

kallisto_all <- 
  read_tsv("./quantifications_ALL/processed_quant.tsv") %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, tpm)

kallisto_ref_imgt <-
  left_join(kallisto_imgt, kallisto_chr, by = c("subject", "locus")) %>%
  left_join(kallisto_all, by = c("subject", "locus")) %>%
  select(subject, locus, imgt = tpm.x, chr = tpm.y, all = tpm)

png("./plots/kallisto_imgt_vs_chr.png", height = 6, width = 10, units = "in", res = 150)
ggplot(kallisto_ref_imgt, aes(imgt, chr)) +
  geom_abline(color = "grey") +
  geom_point(color = "black", alpha = .5) +
  facet_wrap(~locus, scales = "free") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = rel(3.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, margin = margin(t = 7))) +
  labs(x = "kallisto: supplemented index (TPM)", y = "kallisto: REF only (TPM)")
dev.off()

png("./plots/kallisto_imgt_vs_all.png", height = 6, width = 10, units = "in", res = 150)
ggplot(kallisto_ref_imgt, aes(imgt, all)) +
  geom_abline(color = "grey") +
  geom_point(color = "black", alpha = .5) +
  facet_wrap(~locus, scales = "free") +
  ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
                        formula = y ~ x, parse = TRUE, size = rel(3.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, margin = margin(t = 7))) +
  labs(x = "kallisto: supplemented index (TPM)", y = "kallisto: REF + alt haps (TPM)")
dev.off()

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 150)
kallisto_ref_imgt %>%
  gather(index, value, imgt:all) %>%
  ggplot(aes(value, fill = index)) +
  geom_density(alpha = .5) +
  scale_x_continuous(labels = scales::scientific) +
  ggthemes::scale_fill_colorblind() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, margin = margin(t = 10))) +
  facet_wrap(~locus, scales = "free")
dev.off()

# ASE
calc_ase <- function(counts) min(counts)/sum(counts)

pag3f <- mutate(pag, allele = hla_trimnames(allele, 3))

ase_df <- 
  expression_df %>%
  group_by(subject, locus) %>%
  filter(n_distinct(allele) == 2) %>%
  mutate(ase = calc_ase(est_counts)) %>%
  ungroup()
  
ase_error <- 
  calc_genotyping_accuracy(expression_df, pag3f, by_locus = FALSE) %>%
  group_by(subject, locus) %>%
  summarize(error = sum(!correct)) %>%
  ungroup() %>%
  inner_join(ase_df, by = c("subject", "locus"))

png("./plots/ase.png", width = 8, height = 5, units = "in", res = 300)
ggplot(ase_error, aes(factor(error), ase)) +
  geom_jitter(alpha = 1/2, position = position_jitter(width = 0.3)) +
  facet_wrap(~locus) + 
  labs(x = "number of wrong calls in genotype") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        text = element_text(size = rel(4)))
dev.off()

png("./plots/ase_histogram.png", width = 8, height = 5, units = "in", res = 300)
ggplot(ase_df, aes(ase)) +
  geom_density(fill = "grey35", color = NA) +
  facet_wrap(~locus) +
  theme_bw()
dev.off()

# Correlation of expression

## HLA pairs
plot_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
  geom_point(size = .6, alpha = .5) +
  geom_smooth(method = lm) 
}

hla_and_transAct_genes <- gencode_chr_gene %>%
  filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", 
			  "HLA-DRB1", "CIITA"))

pca_expression_files <-
  sprintf("./QTLtools/phenotypes/phenotypes_eur_%d.bed.gz", seq(0, 100, 5)) %>%
  setNames(seq(0, 100, 5))

pca_expression_df <-
  parallel::mclapply(pca_expression_files, 
		     function(i) 
		       read_tsv(i) %>% 
		       inner_join(hla_and_transAct_genes, by = c("id" = "gene_id")) %>%
		       select(gene_name, HG00096:NA20828),
		     mc.cores = 21) %>%
  bind_rows(.id = "PCs") %>%
  gather(subject, value, HG00096:NA20828) %>%
  mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
  spread(gene_name, value) 

phen10 <- filter(pca_expression_df, PCs == 10) %>% select(-PCs)

png("./plots/hlacorrelations.png", width = 10, height = 10, units = "in", res = 300)
pairs_hla_k <-
  ggpairs(select(phen10, -subject, -CIITA), 
          lower = list(continuous = plot_lower), upper = list()) + 
  theme_bw() +
  theme(title = element_text(size = 14))

print(pairs_hla_k, left = .3, bottom = .3)
dev.off()

## Class II transactivator
class_2_trans_df <- 
  phen10 %>%
  select(subject, DQA1, DQB1, DRB1, CIITA) %>%
  gather(locus, value, DQA1, DQB1, DRB1) %>%
  arrange(subject, locus)

png("./plots/trans_activ_corrs.png", width = 10, height = 3.5, units = "in", res = 300)
ggplot(class_2_trans_df, aes(value, CIITA)) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  facet_wrap(~locus) +
  labs(x = NULL) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
               formula = y ~ x, parse = TRUE, size = rel(4))
dev.off()

## Same vs different haplotypes
make_data <- function(locus1, locus2) {
  l1 <- enquo(locus1)
  l2 <- enquo(locus2)

  cis <- select(residuals_by_allele_wide, subject, !!l1, !!l2)
  trans <- cis %>% 
    group_by(subject) %>%
    mutate(!!quo_name(l2) := rev(!!l2)) %>%
    ungroup()
  gene <- select(residuals_gene_level, subject, !!l1, !!l2)
  
  bind_rows(list(Overall = gene, Cis = cis, Trans = trans), .id = "level") %>%
    mutate(level = factor(level, levels = c("Overall", "Cis", "Trans")))
}

make_phase_plot <- function(data, locus1, locus2) {
  ggplot(data, aes_string(locus1, locus2)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = lm, se = FALSE) + 
  facet_wrap(~level, nrow = 1, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(2)),
        strip.text = element_text(size = rel(1.5))) +
  labs(x = paste0("HLA-", locus1), y = paste0("HLA-", locus2)) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
               formula = y ~ x, parse = TRUE, size = rel(5))
}

residuals_by_allele_10pcs <- 
  read_tsv("./phase_hla_alleles/data/hla_allele_expression_10pcs.bed") %>%
  gather(subject, resid, -locus, -hap)

residuals_by_allele_wide <- spread(residuals_by_allele_10pcs, locus, resid) %>%
  select(subject, hap, everything()) %>%
  arrange(subject, hap)

residuals_gene_level <- select(phen10, subject, A, B, C, DQA1, DQB1, DRB1)

png("./plots/a_vs_b.png", height = 3.5, width = 10, units = "in", res = 300)
a_b <- make_data(A, B)
make_phase_plot(a_b, "A", "B")
dev.off()

png("./plots/a_vs_c.png", height = 3.5, width = 10, units = "in", res = 300)
a_c <- make_data(A, C)
make_phase_plot(a_c, "A", "C")
dev.off()

png("./plots/b_vs_c.png", height = 3.5, width = 10, units = "in", res = 300)
b_c <- make_data(B, C)
make_phase_plot(b_c, "B", "C")
dev.off()

png("./plots/dqa_vs_dqb.png", height = 3.5, width = 10, units = "in", res = 300)
dqa_dqb <- make_data(DQA1, DQB1)
make_phase_plot(dqa_dqb, "DQA1", "DQB1")
dev.off()

png("./plots/dqa_vs_drb.png", height = 3.5, width = 10, units = "in", res = 300)
dqa_drb <- make_data(DQA1, DRB1)
make_phase_plot(dqa_drb, "DQA1", "DRB1")
dev.off()

# PCA vs PEER correction
cors_pca <- 
  pca_expression_df %>%
  group_by(PCs) %>%
  summarize(AxB = cor(A, B), 
            AxC = cor(A, C), 
	    BxC = cor(B, C),
	    DQA1xDQB1 = cor(DQA1, DQB1), 
	    DQA1xDRB1 = cor(DQA1, DRB1), 
	    DQB1xDRB1 = cor(DQB1, DRB1),  
	    DQA1xCIITA = cor(DQA1, CIITA),
	    DQB1xCIITA = cor(DQB1, CIITA), 
	    DRB1xCIITA = cor(DRB1, CIITA)) %>%
  gather(gene_pair, correlation, -1)

peer_expression_files <- 
  list.files("./QTLtools/peer/phenotypes", pattern = "\\.bed\\.gz$", full.names = TRUE)
 
names(peer_expression_files) <- sub("^.+eur_(\\d+)\\.bed\\.gz$", "\\1", peer_expression_files)

peer_expression_df <-
  parallel::mclapply(peer_expression_files, 
		     function(i) 
		       read_tsv(i) %>% 
		       inner_join(hla_and_transAct_genes, by = c("id" = "gene_id")) %>%
		       select(gene_name, HG00096:NA20828),
		     mc.cores = length(peer_expression_files)) %>%
  bind_rows(.id = "K") %>%
  gather(subject, value, HG00096:NA20828) %>%
  mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
  spread(gene_name, value) 

cors_peer <- 
  peer_expression_df %>%
  group_by(K) %>%
  summarize(AxB = cor(A, B), 
            AxC = cor(A, C), 
	    BxC = cor(B, C),
	    DQA1xDQB1 = cor(DQA1, DQB1), 
	    DQA1xDRB1 = cor(DQA1, DRB1), 
	    DQB1xDRB1 = cor(DQB1, DRB1),  
	    DQA1xCIITA = cor(DQA1, CIITA),
	    DQB1xCIITA = cor(DQB1, CIITA), 
	    DRB1xCIITA = cor(DRB1, CIITA)) %>%
  gather(gene_pair, correlation, -1)

cors_data <- 
  left_join(cors_pca, cors_peer, by = c("gene_pair", "PCs" = "K")) %>%
  rename(covariates = PCs, PCA = correlation.x, PEER = correlation.y) %>%
  gather(method, correlation, PCA:PEER) %>%
  mutate(covariates = as.integer(covariates),
         gene_pair = factor(gene_pair, levels = c("AxB", "AxC", "BxC",
                                                  "DQA1xDQB1", "DQA1xDRB1", "DQB1xDRB1",
                                                  "DQA1xCIITA", "DQB1xCIITA", "DRB1xCIITA"))) %>%
  arrange(covariates, gene_pair, method)

png("./plots/correlation_decrease.png", width = 10, height = 5, units = "in", res = 300)
ggplot(cors_data, aes(covariates, correlation, color = method)) +
  geom_point() +
  ggsci::scale_color_npg() +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  facet_wrap(~gene_pair) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
  labs(x = "Number of PCs/factors")
dev.off()