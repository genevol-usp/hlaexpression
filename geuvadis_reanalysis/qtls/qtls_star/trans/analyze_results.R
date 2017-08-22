devtools::load_all("~/hlaseqlib")
library(tidyverse)

gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
  select(gene_name, gene_id)

hits <- 
  read_delim("./fdr_pvals.txt", delim = " ", col_types = "ciiciid--",
	     col_names = c("gene_id", "chr_gene", "start_gene", "variant_id", 
			   "variant_chr", "variant_pos", "nom_pval")) %>%
  inner_join(gencode_hla, by = "gene_id") %>%
  mutate(log10_nom_pval = -log10(nom_pval)) %>%
  arrange(desc(log10_nom_pval)) %>%
  select(gene_name, variant_id, variant_chr, variant_pos, log10_nom_pval) 

write_tsv(hits, "./trans_results.tsv")
