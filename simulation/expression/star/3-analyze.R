devtools::load_all("~/hlatools")
devtools::load_all("~/genomicRutils")
library(tidyverse)

gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_id, gene_name)

genos <- read_tsv("../../../simulation/ground_truth_files/genos.tsv")

samples <- sprintf("sample_%02d", 1:50)

quant <- 
  file.path("./quantifications_1", samples, "quant.sf") %>%
  setNames(samples) %>%
  map_df(~read_tsv(., col_types = "c--dd") %>% 
	 filter(grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", Name)) %>%
	 mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
	 select(locus, allele = Name, est_counts = NumReads, tpm = TPM),
       .id = "subject")

quant_typed <- 
  left_join(quant, gencode_hla, by = c("locus" = "gene_name")) %>%
  hla_genotype_dt(th = 0.05) %>%
  mutate(allele = hla_trimnames(gsub("IMGT_|_s\\d", "", allele), 3)) %>%
  select(subject, locus, allele)

calc_genotyping_accuracy(quant_typed, genos)

typings <- calc_genotyping_accuracy(quant_typed, genos, by_locus = FALSE)
  
