devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(Biostrings)
library(tidyverse)

genome <- readDNAStringSet("/home/vitor/gencode_data/GRCh38.primary_assembly.genome.fa.gz")
chr6 <- genome[names(genome) == "chr6 6"][[1]]
mhc <- chr6[25e6:35e6] 

gen <- readDNAStringSet("/home/vitor/IMGTHLA/hla_gen.fasta")
names(gen) <- sub("^[^ ]+ ([^ ]+).*$", "\\1", names(gen))

gen_plus <- gencode_hla %>%
    filter(strand == "+") %>%
    pull(gene_name) %>%
    sub("HLA-", "", .) %>%
    paste(collapse = "|") %>%
    paste0("^(", ., ")\\*") %>%
    grep(names(gen)) %>%
    gen[.]
    
gen_minus <- gencode_hla %>%
    filter(strand == "-") %>%
    pull(gene_name) %>%
    sub("HLA-", "", .) %>%
    paste(collapse = "|") %>%
    paste0("^(", ., ")\\*") %>%
    grep(names(gen)) %>%
    gen[.] %>%
    reverseComplement()
 
hlagen <- c(gen_plus, gen_minus)

ref_alleles <- matchPDict(hlagen, chr6) %>%
    as.list() %>%
    map_df(as.data.frame, .id = 'ref_allele') %>%
    mutate(locus = sub("^([^*]+).+$", "HLA-\\1", ref_allele)) %>%
    arrange(start) %>%
    select(locus, ref_allele) %>%
    as_tibble()

write_tsv(ref_alleles, "./ref_alleles.tsv")
