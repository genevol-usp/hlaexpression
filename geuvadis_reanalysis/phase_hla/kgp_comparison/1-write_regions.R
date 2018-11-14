devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

hla_genes <- c(gencode_hla$gene_name, "HLA-E", "HLA-G") 

hla_annot <- filter(gencode_chr_gene, gene_name %in% hla_genes)

hla_annot %>%
    pull(gene_name) %>%
    sub("HLA-", "", .) %>%
    writeLines("./imgt_loci.txt")

bed <- select(hla_annot, chr, start, end)

write_tsv(bed, "./hla_regions.bed", col_names = FALSE)
