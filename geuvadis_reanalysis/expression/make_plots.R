devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
library(ggpmisc)

# functions
calc_ase <- function(counts) min(counts)/sum(counts)

plot_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
  geom_point(size = .6, alpha = .5) +
  geom_smooth(method = lm) 
}

scatter_plot <- function(df, x_var, y_var) {
  ggplot(df, aes_string(x_var, y_var)) +
    geom_abline() +
    geom_point() +
    facet_wrap(~locus, scales = "free") +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                          formula = y ~ x, parse = TRUE, size = 6) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))
}

scatter_plot_color <- function(df, x_var, y_var, alpha_var) {
  ggplot(df, aes_string(x_var, y_var, alpha = alpha_var)) +
    geom_abline() +
    geom_point() +
    facet_wrap(~locus, scales = "free") +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                          formula = y ~ x, parse = TRUE, size = 6) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16)) +
    labs(alpha = "distance to ref")
}

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
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) + 
  labs(x = paste0("HLA-", locus1), y = paste0("HLA-", locus2)) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
               formula = y ~ x, parse = TRUE, size = 6)
}

# data
geuvadis_ids <- geuvadis_info %>%
  filter(kgp_phase3 == 1, pop != "YRI") %>%
  select(subject = ena_id, name)

library_size <- 
  read_delim("../data/sample_info/library_size.txt", delim = " ") %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject, total = n) %>%
  arrange(subject)

kallisto_aligned <- read_tsv("./kallisto/aligned_reads.tsv")
star_aligned <- read_tsv("./star/aligned_reads.tsv")

reads_df <- 
  left_join(kallisto_aligned, star_aligned, by = "subject") %>%
  rename(kallisto = V1.x, STAR_salmon = V1.y) %>%
  left_join(library_size, by = "subject") %>%
  gather("source", "n_reads", -1) %>%
  mutate(source = factor(source, levels = c("total", "kallisto", "STAR_salmon"))) %>%
  arrange(subject, source)

gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")))

gencode12 <-
  "/home/vitor/gencode_data/gencode.v12.annotation.gtf.gz" %>%
  get_gencode_coords(feature = "gene") %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_name, gene_id)

gencode19 <-
  "/home/vitor/gencode_data/gencode.v19.annotation.gtf.gz" %>%
  get_gencode_coords(feature = "gene") %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_name, gene_id)

distances <- read_tsv("../../simulation/data/distances_to_reference.tsv")

hlatx <- 
  hla_readtx("./hlaTX/results/results_TX.txt") %>%
  mutate(locus = sub("^([^_]+).*$", "HLA-\\1", allele)) %>%
  inner_join(geuvadis_ids) %>%
  select(subject = name, locus, est_counts) %>%
  group_by(subject, locus) %>%
  summarize(counts = sum(est_counts)) %>%
  ungroup()

kallisto <- 
  read_tsv("./kallisto/quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1")) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, allele, est_counts, tpm) %>%
  mutate(locus = paste0("HLA-", locus),
	 allele = sub("IMGT_", "", allele)) %>%
  left_join(distances, by = c("locus", "allele")) %>%
  arrange(subject, locus, allele)

kallisto_gene <- kallisto %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
  ungroup()

kallisto_10pc <-
  "../qtls/qtls_kallisto/qtltools_correction/phenotypes/phenotypes_eur_10.bed.gz" %>%
  read_tsv() %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

star_10pc <-
  "../qtls/qtls_star/phenotypes/phenotypes_eur_10.bed.gz" %>%
  read_tsv() %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

star_gene <- 
  read_tsv("./star/quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, allele, est_counts, tpm) %>%
  mutate(allele = sub("IMGT_", "", allele)) %>%
  left_join(distances, by = c("locus", "allele")) %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
  ungroup()

geuvadis <- 
  read_tsv("../data/quantifications/peer/published/phenotypes_eur_10.bed.gz") %>%
  inner_join(gencode12, by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

geuvadis_new <-
  read_tsv("../data/quantifications/peer/new/phenotypes_eur_10.bed.gz") %>%
  inner_join(gencode19, by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

kallisto_fpkms_10pcs <- 
  read_tsv("../qtls/qtls_kallisto/peer_correction/phenotypes_fpkm/phenotypes_eur_10.bed.gz") %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

star_fpkms_10pcs <-
  read_tsv("../qtls/qtls_star/phenotypes_fpkm_peer/phenotypes_eur_10.bed.gz") %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

kallisto_chr <- 
  "../qtls/qtls_kallisto/qtltools_correction/phenotypes_chr/phenotypes_eur_10.bed.gz" %>%
  read_tsv() %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

star_chr <- 
  "../qtls/qtls_star/phenotypes_chr/phenotypes_eur_10.bed.gz" %>%
  read_tsv() %>%
  inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
  select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
  gather(subject, resid, -locus)

quant_data <- 
  left_join(kallisto_gene, hlatx, by = c("subject", "locus")) %>%
  select(subject, locus, dist.kallisto = dist, 
	 est_counts.kallisto.imgt = est_counts, 
	 est_counts.hlatx.imgt = counts, tpm.kallisto.imgt = tpm) %>%
  left_join(library_size, by = "subject") %>%
  mutate(est_counts.kallisto.imgt.byLibSize = est_counts.kallisto.imgt/total * 1e6, 
	 est_counts.hlatx.imgt.byLibSize = est_counts.hlatx.imgt/total * 1e6) %>%
  select(-total) %>%
  left_join(star_gene, by = c("subject", "locus")) %>%
  rename(est_counts.star.imgt = est_counts, tpm.star.imgt = tpm, dist.star = dist) %>%
  left_join(kallisto_10pc, by = c("subject", "locus")) %>%
  rename(resid.kallisto.imgt.tpm = resid) %>%
  left_join(star_10pc, by = c("subject", "locus")) %>%
  rename(resid.star.imgt.tpm = resid) %>%
  left_join(kallisto_chr, by = c("subject", "locus")) %>%
  rename(resid.kallisto.chr.tpm = resid) %>%
  left_join(kallisto_fpkms_10pcs, by = c("subject", "locus")) %>%
  rename(resid.kallisto.imgt.fpkm = resid) %>%
  left_join(star_fpkms_10pcs, by = c("subject", "locus")) %>%
  rename(resid.star.imgt.fpkm = resid) %>%
  left_join(star_chr, by = c("subject", "locus")) %>%
  rename(resid.star.chr.tpm = resid) %>%
  left_join(geuvadis, by = c("subject", "locus")) %>%
  rename(resid.geuvadis.old = resid) %>%
  left_join(geuvadis_new, by = c("subject", "locus")) %>%
  rename(resid.geuvadis.new = resid)

kallisto_rsq <- 
  quant_data %>%
  group_by(locus) %>%
  mutate(qt_dist = ntile(dist.kallisto, 4)) %>%
  filter(qt_dist %in% c(1, 4)) %>%
  group_by(locus, qt_dist) %>%
  summarize(rsq = summary(lm(resid.kallisto.imgt.tpm~resid.kallisto.chr.tpm))$adj.r.square) %>%
  ungroup()

star_rsq <- 
  quant_data %>%
  group_by(locus) %>%
  mutate(qt_dist = ntile(dist.star, 4)) %>%
  filter(qt_dist %in% c(1, 4)) %>%
  group_by(locus, qt_dist) %>%
  summarize(rsq = summary(lm(resid.star.imgt.tpm~resid.star.chr.tpm))$adj.r.square) %>%
  ungroup()

kallisto_imgt_tpm <- select(kallisto_gene, -est_counts, -dist)

kallisto_chr_tpm <- 
  read_tsv("./kallisto/quantifications_CHR/processed_quant.tsv") %>%
  mutate(locus = paste0("HLA-", locus)) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, tpm) %>%
  arrange(subject, locus)

kallisto_all_tpm <- 
  read_tsv("./kallisto/quantifications_ALL/processed_quant.tsv") %>%
  mutate(locus = paste0("HLA-", locus)) %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, tpm) %>%
  arrange(subject, locus)

kallisto_tpm_df <- 
  list(imgt = kallisto_imgt_tpm, chr = kallisto_chr_tpm, all = kallisto_all_tpm) %>%
  bind_rows(.id = "index")

kallisto_tpm_rsq <- 
  select(kallisto_gene, -est_counts) %>%
  left_join(kallisto_chr_tpm, by = c("subject", "locus"), 
	    suffix = c(".kallisto.imgt", ".kallisto.chr")) %>%
  group_by(locus) %>%
  mutate(qt_dist = ntile(dist, 4)) %>%
  filter(qt_dist %in% c(1, 4)) %>%
  group_by(locus, qt_dist) %>%
  summarize(rsq = summary(lm(tpm.kallisto.imgt~tpm.kallisto.chr))$adj.r.square) %>%
  ungroup()

star_imgt_tpm <- star_gene %>%
  select(subject, locus, dist, tpm)

star_chr_tpm <- 
  read_tsv("./star/quantifications_CHR/processed_quant.tsv") %>%
  inner_join(geuvadis_ids, by = "subject") %>%
  select(subject = name, locus, tpm)

star_tpm_rsq <- 
  left_join(star_imgt_tpm, star_chr_tpm, by = c("subject", "locus"), 
	    suffix = c(".kallisto.imgt", ".kallisto.chr")) %>%
  group_by(locus) %>%
  mutate(qt_dist = ntile(dist, 4)) %>%
  filter(qt_dist == 1L | qt_dist == 4L) %>%
  group_by(locus, qt_dist) %>%
  summarize(rsq = summary(lm(tpm.kallisto.imgt~tpm.kallisto.chr))$adj.r.square) %>%
  ungroup()

rsq_df <-
  left_join(kallisto_rsq, star_rsq, by = c("locus", "qt_dist"), 
	    suffix = c(".kallisto.pca", ".star.pca")) %>%
  left_join(kallisto_tpm_rsq, by = c("locus", "qt_dist")) %>%
  rename(rsq.kallisto.tpm = rsq) %>%
  left_join(star_tpm_rsq, by = c("locus", "qt_dist")) %>%
  rename(rsq.star.tpm = rsq)

write_tsv(rsq_df, "./plots/rsq_by_quartile.tsv")

pag3f <- mutate(pag, allele = hla_trimnames(allele, 3))

ase_df <- 
  kallisto %>%
  group_by(subject, locus) %>%
  filter(n_distinct(allele) == 2) %>%
  summarize(ase = calc_ase(est_counts)) %>%
  ungroup()
  
ase_error <- 
  calc_genotyping_accuracy(kallisto, pag3f, by_locus = FALSE) %>%
  group_by(subject, locus) %>%
  summarize(error = sum(!correct)) %>%
  ungroup() %>%
  inner_join(ase_df, by = c("subject", "locus"))

hla_and_transAct_genes <- gencode_chr_gene %>%
  filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", 
			  "HLA-DRB1", "CIITA"))

pcs <- c(seq(0, 30, 5), seq(40, 100, 10))

pca_kallisto_df <-
  sprintf("../qtls/qtls_kallisto/qtltools_correction/phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
  setNames(pcs) %>%
  parallel::mclapply(function(i) 
		       read_tsv(i, progress = FALSE) %>% 
		       inner_join(hla_and_transAct_genes, by = c("id" = "gene_id")) %>%
		       select(gene_name, HG00096:NA20828),
		     mc.cores = length(pcs)) %>%
  bind_rows(.id = "PCs") %>%
  gather(subject, value, HG00096:NA20828) %>%
  mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
  spread(gene_name, value) 

pca_star_df <-
  sprintf("../qtls/qtls_star/phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
  setNames(pcs) %>%
  parallel::mclapply(function(i) 
		       read_tsv(i, progress = FALSE) %>% 
		       inner_join(hla_and_transAct_genes, by = c("id" = "gene_id")) %>%
		       select(gene_name, HG00096:NA20828),
		     mc.cores = length(pcs)) %>%
  bind_rows(.id = "PCs") %>%
  gather(subject, value, HG00096:NA20828) %>%
  mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
  spread(gene_name, value) 

phen10 <- filter(pca_kallisto_df, PCs == 10) %>% select(-PCs)

class_2_trans_df <- 
  phen10 %>%
  select(subject, DQA1, DQB1, DRB1, CIITA) %>%
  gather(locus, value, DQA1, DQB1, DRB1) %>%
  arrange(subject, locus)

residuals_by_allele_10pcs <- 
  read_tsv("./kallisto/phase_hla_alleles/data/hla_allele_expression_10pcs.bed") %>%
  gather(subject, resid, -locus, -hap)

residuals_by_allele_wide <- 
  spread(residuals_by_allele_10pcs, locus, resid) %>%
  select(subject, hap, everything()) %>%
  arrange(subject, hap)

residuals_gene_level <- select(phen10, subject, A, B, C, DQA1, DQB1, DRB1)

cors_pca_kallisto <- 
  pca_kallisto_df %>%
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

cors_pca_star <- 
  pca_star_df %>%
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
  list.files("../qtls/qtls_kallisto/peer_correction/phenotypes", pattern = "\\.bed\\.gz$", full.names = TRUE)
 
names(peer_expression_files) <- sub("^.+eur_(\\d+)\\.bed\\.gz$", "\\1", peer_expression_files)

peer_expression_df <-
  parallel::mclapply(peer_expression_files, 
		     function(i) 
		       read_tsv(i, progress = FALSE) %>% 
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
  left_join(cors_pca_kallisto, cors_pca_star, by = c("PCs", "gene_pair")) %>%
  rename(kallisto.PCA = correlation.x, star_salmon.PCA = correlation.y) %>%
  left_join(cors_peer, by = c("gene_pair", "PCs" = "K")) %>%
  rename(covariates = PCs, kallisto.PEER = correlation) %>%
  gather(method, correlation, -covariates, -gene_pair) %>%
  mutate(covariates = as.integer(covariates),
         gene_pair = factor(gene_pair, levels = c("AxB", "AxC", "BxC",
                                                  "DQA1xDQB1", "DQA1xDRB1", "DQB1xDRB1",
                                                  "DQA1xCIITA", "DQB1xCIITA", "DRB1xCIITA"))) %>%
  arrange(covariates, gene_pair, method)

# plots
png("./plots/library_sizes.png", height = 4, width = 10, units = "in", res = 150)
ggplot(reads_df) +
  geom_line(aes(x = reorder(subject, n_reads, FUN = "max"), y = n_reads, 
                color = source, group = source), size = 1.1) +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()

png("./plots/kallisto_vs_hlatx.png", height = 6, width = 10, units = "in", res = 150)
scatter_plot(quant_data, "est_counts.kallisto.imgt.byLibSize", "est_counts.hlatx.imgt.byLibSize") +
  labs(x = "Counts per million reads (kallisto)", 
       y = "Counts per million reads (GEM-based method)")
dev.off()

png("./plots/kallisto_vs_star.png", width = 10, height = 6, units = "in", res = 300)
scatter_plot(quant_data, "resid.kallisto.imgt.tpm", "resid.star.imgt.tpm") +
  labs(x = "PCA-corrected TPM (kallisto)", 
       y = "PCA-corrected TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_vs_geuvadis.png", width = 10, height = 6, units = "in", res = 300)
scatter_plot_color(quant_data, "resid.kallisto.imgt.fpkm", "resid.geuvadis.old", "dist.kallisto") +
  labs(x = "PEER-corrected FPKM (kallisto)", y = "PEER-corrected FPKM (GEUVADIS)")
dev.off()

png("./plots/star_vs_geuvadis.png", width = 10, height = 6, units = "in", res = 300)
scatter_plot_color(quant_data, "resid.star.imgt.fpkm", "resid.geuvadis.old", "dist.star") +
  labs(x = "PEER-corrected FPKM (STAR)", y = "PEER-corrected FPKM (GEUVADIS)")
dev.off()

png("./plots/kallisto_vs_geuvadis_new.png", width = 10, height = 6, units = "in", res = 300)
scatter_plot_color(quant_data, "resid.kallisto.imgt.fpkm", "resid.geuvadis.new", "dist.kallisto") +
  labs(x = "PEER-corrected FPKM (kallisto)", y = "PEER-corrected FPKM (GEUVADIS)")
dev.off()

png("./plots/star_vs_geuvadis_new.png", width = 10, height = 6, units = "in", res = 300)
scatter_plot_color(quant_data, "resid.star.imgt.fpkm", "resid.geuvadis.new", "dist.star") +
  labs(x = "PEER-corrected FPKM (STAR)", y = "PEER-corrected FPKM (GEUVADIS)")
dev.off()

png("./plots/kallisto_imgt_vs_chr.png", height = 6, width = 10, units = "in", res = 150)
scatter_plot_color(quant_data, "resid.kallisto.imgt.tpm", "resid.kallisto.chr.tpm", "dist.kallisto") +
  labs(x = "PCA-corrected TPM (kallisto-IMGT)", y = "PCA-corrected (kallisto REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_chr.png", height = 6, width = 10, units = "in", res = 150)
scatter_plot_color(quant_data, "resid.star.imgt.tpm", "resid.star.chr.tpm", "dist.star") +
  labs(x = "PCA-corrected TPM (STAR-IMGT)", y = "PCA-corrected (STAR REF chromosomes)")
dev.off()

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 150)
ggplot(kallisto_tpm_df, aes(tpm, fill = index)) +
  geom_density(alpha = 1/2) +
  scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
  ggthemes::scale_fill_colorblind() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, margin = margin(t = 10))) +
  facet_wrap(~locus, scales = "free")
dev.off()

png("./plots/ase.png", width = 8, height = 5, units = "in", res = 300)
ggplot(ase_error, aes(factor(error), ase)) +
  ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = .75, alpha = 1/2) +
  facet_wrap(~locus) + 
  labs(x = "number of wrong calls in genotype") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
dev.off()

png("./plots/ase_histogram.png", width = 8, height = 4, units = "in", res = 300)
ggplot(ase_df, aes(ase)) +
  geom_density(fill = "grey35", color = NA) +
  facet_wrap(~locus) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
dev.off()

png("./plots/hlacorrelations.png", width = 8, height = 8, units = "in", res = 300)
pairs_hla_k <-
  ggpairs(select(phen10, -subject, -CIITA), 
          lower = list(continuous = plot_lower), upper = list()) + 
  theme_bw() +
  theme(title = element_text(size = 14))

print(pairs_hla_k, left = .3, bottom = .3)
dev.off()

png("./plots/trans_activ_corrs.png", width = 10, height = 3.5, units = "in", res = 300)
ggplot(class_2_trans_df, aes(value, CIITA)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) + 
  facet_wrap(~locus) +
  labs(x = NULL) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
               formula = y ~ x, parse = TRUE, size = 6)
dev.off()

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

png("./plots/correlation_decrease.png", width = 10, height = 5, units = "in", res = 300)
ggplot(cors_data, aes(covariates, correlation, color = method)) +
  geom_point() +
  ggsci::scale_color_npg() +
  scale_x_continuous(breaks = pcs) +
  facet_wrap(~gene_pair) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
  labs(x = "Number of PCs/factors")
dev.off()
