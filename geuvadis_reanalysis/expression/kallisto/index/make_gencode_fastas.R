library(Biostrings)
devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

imgt <- readDNAStringSet("imgt_index.fa") 

gencode <- readDNAStringSet("/home/vitor/gencode_data/gencode.v25.transcripts.fa.gz")
  
names(gencode) <- sub("^([^|]+).*$", "\\1", names(gencode))  

loci_in_index <- unique(imgt_to_gname(names(imgt)))

gencode_chr_imgt <-
  gencode_chr_tx %>%
  filter(gene_name %in% loci_in_index) %>%
  arrange(gene_name, tx_id)

gencode_no_imgt <- gencode[! names(gencode) %in% gencode_chr_imgt$tx_id] 
  
gencode_imgt <- c(gencode_no_imgt, imgt) 

writeXStringSet(gencode, "./gencode.v25.CHR.transcripts.fa")
writeXStringSet(gencode_imgt, "./gencode.v25.CHR.IMGT.transcripts.fa")
writeXStringSet(gencode_no_imgt, "./gencode.v25.CHR.transcripts.noIMGT.fa")
