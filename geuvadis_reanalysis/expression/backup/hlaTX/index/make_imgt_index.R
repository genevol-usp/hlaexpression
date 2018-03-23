devtools::load_all("/home/vitor/hlatools")

make_index <- function(seqs_list)
  seqs_list %>%
  purrr::flatten() %>%
  purrr::set_names(names(.) %>% 
		   allele_to_group() %>%
		   hla_trimnames(3) %>% 
		   stringr::str_replace_all("\\*|:", "_") %>%
		   stringr::str_c("HLA_", .)) %>%
  .[! duplicated(names(.))]

alignments <- "/home/vitor/imgt_data/alignments"

df <- 
  tibble::tribble(
    ~locus,  ~ars_exons,
    "A",     2:3,
    "B",     2:3,
    "C",     2:3,
    "DQA1",  2L,
    "DQB1",  2L,
    "DRB",   2L)

whole_seqs <- vector("list", nrow(df))
ars_seqs <- vector("list", nrow(df))

for (i in seq_len(nrow(df))) {

  locus <- df$locus[[i]]
  message(stringr::str_c("processing locus ", locus))

  ars_seqs[[i]] <- 
    hla_parsenuc(locus, alignments, exons = df$ars_exons[[i]]) %>%
    .[stringr::str_detect(rownames(.), "^[A-Z]{1,3}[^2-9]*\\*.+[^N]$"), ] %>%
    hla_list(rm_asterisks = TRUE)

  whole_seqs[[i]] <-  
    hla_parsenuc(locus, alignments) %>%
    .[stringr::str_detect(rownames(.), "^[A-Z]{1,3}[^2-9]*\\*.+[^N]$"), ] %>%
    hla_infer(cores = 40) %>%
    hla_list(rm_asterisks = TRUE)
}
  
short_index <- make_index(ars_seqs)

long_index <- make_index(whole_seqs)

seqinr::write.fasta(short_index, names(short_index), "index_short.fasta")	      
seqinr::write.fasta(long_index, names(long_index), "index_all_exons.fasta")	    

system("./run_gem_indexer.sh")
