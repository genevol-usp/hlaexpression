devtools::load_all("~/hlaseqlib")
library(tidyverse)

read_gencode <- function(gtf) {

    read_tsv(gtf, comment = "##", col_names = FALSE,
	     col_types = "c-cii-c-c", progress = FALSE) %>%
    select(chr = X1, feature = X3, start = X4, end = X5, strand = X7, X9)
}

gencode25 <- read_gencode("~/gencode_data/gencode.v25.annotation.gtf.gz") %>%
    filter(chr == "chr6", strand %in% c("+", "-"), feature != "transcript")

qtls <- "../3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(gene_name, rank, var_id, var_chr, var_from) %>%
    mutate(var_chr = paste0("chr", var_chr))

x <- 
    full_join(qtls, gencode25, by = c("var_chr" = "chr")) %>%
    distinct(gene_name, var_id, feature, start, end, .keep_all = TRUE) %>%
    mutate(b = var_from >= start & var_from <= end,
	   ups_d = ifelse(strand == "+", var_from - start, end - var_from),
	   dow_d = ifelse(strand == "+", var_from - end, start - var_from),
	   gene_ups = ups_d, gene_dow = dow_d,
	   ups_d = ifelse(b, 0, ups_d),
	   dow_d = ifelse(b, 0, dow_d)) %>%
    select(-b, -var_chr) %>%
    gather(ori, distance, ups_d, dow_d) %>% 
    group_by(gene_name, rank, ori) %>%
    filter(abs(distance) == min(abs(distance))) %>%
    ungroup() %>%  
    mutate(X9 = gsub("\"", "", X9),
	   gene = sub("^.+gene_name ([^;]+).*$", "\\1", X9)) %>%
    select(-X9) %>%
    group_by(gene_name, rank, gene) %>%
    mutate(exon = as.integer(all(c("gene", "exon") %in% feature)),
	   utr = as.integer(all(c("exon", "UTR") %in% feature))) %>%
    filter(distance != 0 | (distance == 0 & ori == "ups_d")) %>%
    ungroup() %>% 
    arrange(gene_name, rank, ori)

x %>% 
    print(n = Inf)

