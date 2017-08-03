library(tidyverse)

allele_dist <- read_tsv("./data/distances_to_reference.tsv") %>%
  mutate(locus = sub("^HLA-", "", locus))

samples <- tibble(subject = readLines("./data/samples.txt"),
                  code = sprintf("sample_%02d", 1:50))

index <- Biostrings::readDNAStringSet("./data/polyester_index.fa")

ground_truth <- 
  read_tsv("./data/phenotypes.tsv") %>%
  mutate(target_id = names(index)) %>%
  filter(grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", target_id)) %>%
  gather(subject, true_counts, -target_id) %>%
  inner_join(samples, by = "subject") %>%
  select(subject = code, target_id, true_counts) %>%
  mutate(locus = sub("^IMGT_([^*]+).+$", "\\1", target_id)) %>%
  group_by(subject, locus) %>%
  summarize(true_counts = sum(true_counts))

quants_imgt <- 
  read_tsv("./expression/kallisto/quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1")) %>%
  mutate(allele = sub("IMGT_", "", allele)) %>%
  left_join(allele_dist, by = c("locus", "allele")) %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts), dist = mean(dist))

quants_chr <- 
  read_tsv("./expression/kallisto/quantifications_CHR/processed_quant.tsv") %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts))

quants_all <- 
  read_tsv("./expression/kallisto/quantifications_ALL/processed_quant.tsv") %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts))

quant_data <-
  left_join(quants_chr, quants_all, by = c("subject", "locus"),
            suffix = c("_chr", "_all")) %>%
  left_join(quants_imgt, by = c("subject", "locus")) %>%
  rename(est_counts_imgt = est_counts) %>%
  gather(index, counts, est_counts_chr:est_counts_imgt) %>%
  left_join(ground_truth, by = c("subject", "locus")) %>%
  mutate(index = sub("est_counts_", "", index), 
         prop_mapped = counts/true_counts)

quant_star <- 
  read_tsv("./expression/star/quantifications_2/processed_quant.tsv") %>%
  filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  group_by(subject, locus) %>%
  summarize(est_counts = sum(est_counts)) %>%
  mutate(locus = sub("HLA-", "", locus)) %>%
  left_join(ground_truth, by = c("subject", "locus"))

png("./plots/hlasimul.png", width = 12, height = 6, units = "in", res = 300)
ggplot(quant_data, aes(dist, prop_mapped, color = index)) +
  geom_point() +
  geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
            alpha = 0.4, size = 1.5) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(breaks = seq(0, 1.4, 0.4)) +
  facet_wrap(~locus, scales = "free_x") +
  ggsci::scale_color_aaas(labels = c(chr = "Ref chromosomes",
                                    all = stringr::str_wrap("Ref chromosomes + Alternate haplotypes", 20),
                                    imgt = stringr::str_wrap("Ref chromosomes + HLA diversity", 20))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "proportion of sites with mismatches to the reference allele", 
       y = "proportion of reads recovered")
dev.off()

png("./plots/star_simul.png", width = 8, height = 5, units = "in", res = 300)
ggplot(quant_star, aes(est_counts, true_counts)) +
  geom_abline() +
  geom_point(alpha = .5) +
  facet_wrap(~locus, scales = "free") +
  theme_bw() +
  labs(x = "estimated counts", y = "true counts")
dev.off()
