library(tidyverse)

overlap_elements <- function(pos, element_df) {

    element_df %>%
	mutate(i = map2_lgl(X2, X3, ~between(pos, .x, .y))) %>%
	filter(i) %>%
	pull(X4) %>%
	unique() %>%
	paste(collapse = "/")
}

dnase <- read_tsv("./DNase.ENCODE.chr6.hg38.bed", col_names = FALSE) %>%
    select(X1:X3) %>%
    mutate(X4 = "DHS") %>%
    distinct()

tf <- read_tsv("./TF.ENCODE.chr6.hg38.bed", col_names = FALSE) %>%
    distinct()

seg <- read_tsv("./Seg.ENCODE.chr6.hg38.bed", col_names = FALSE) %>%
    select(X1:X4) %>%
    distinct()

hm <- read_tsv("./HM.ENCODE.chr6.hg38.bed", col_names = FALSE) %>%
    distinct()

qtl <- 
    read_tsv("./hla.qtls.bed", col_names = FALSE) %>%
    mutate(tf = map_chr(X2, overlap_elements, tf),
	   dhs = map_chr(X2, overlap_elements, dnase),
	   chrom_state = map_chr(X2, overlap_elements, seg),
	   histone_marks = map_chr(X2, overlap_elements, hm)) %>%
    separate(X4, c("locus", "rsid", "rank"), sep = ":", convert = TRUE) %>%
    select(locus, rank, rsid, pos = X3, tf, dhs, chrom_state, histone_marks) %>%
    arrange(locus, rank) %>%
    mutate_at(vars(tf, dhs, chrom_state, histone_marks), ~ifelse(. == "", NA, .))

write_tsv(qtl, "./hla.qtl.functional.elements.tsv")
