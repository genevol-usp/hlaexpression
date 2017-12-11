devtools::load_all("~/hlaseqlib")
library(data.table)

phenotypes <- 
  fread("zcat < ../../qtls/qtls_kallisto/qtltools_correction/phenotypes/phenotypes_eur_10.bed.gz")

setDT(gencode_chr_gene)

gencode_hla <- 
  gencode_chr_gene[gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))
		 ][, `:=`(start = start - 5e6, end = end + 5e6)]
setkey(gencode_hla, chr, start, end)

cis5m_genes <- 
  foverlaps(gencode_chr_gene, gencode_hla, type = "within", nomatch = 0L
	  )[gene_id != i.gene_id, .(gene_id, gene_name, i.gene_id, i.gene_name)]

