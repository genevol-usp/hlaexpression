devtools::load_all("/home/vitor/genomicRutils")

imgt <- Biostrings::readDNAStringSet("imgt_index.fa") 

gencode <- 
  Biostrings::readDNAStringSet("~/gencode_data/gencode.v26.transcripts.fa.gz")
  
names(gencode) <- sub("^([^|]+).*$", "\\1", names(gencode))  

loci_in_index <- unique(imgt_to_gname(names(imgt)))

gencode_chr_imgt <-
  gencode_chr_tx %>%
  dplyr::filter(gene_name %in% loci_in_index) %>%
  dplyr::arrange(gene_name, tx_id)

gencode_no_imgt <- gencode[! names(gencode) %in% gencode_chr_imgt$tx_id] 
  
gencode_imgt <- c(gencode_no_imgt, imgt) 

Biostrings::writeXStringSet(gencode, "./gencode.v26.CHR.transcripts.fa")
Biostrings::writeXStringSet(gencode_imgt, "./gencode.v26.CHR.IMGT.transcripts.fa")
Biostrings::writeXStringSet(gencode_no_imgt, "./gencode.v26.CHR.transcripts.noIMGT.fa")
