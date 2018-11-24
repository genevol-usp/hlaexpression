devtools::load_all("~/Libraries/hlaseqlib")
library(tidyverse)

phase_out <- readLines("./phase.out")

phase_best <- grep("(BEGIN|END) BESTPAIRS1", phase_out) %>% 
    {(.[[1]] + 1):(.[[2]] - 1)} %>%
    phase_out[.] 
    
best_list <- split(phase_best, grepl("#", phase_best) %>% cumsum)

loci_df <- read_tsv("./input_loci.tsv") %>%
    select(locus, pos)

hla_codes <- read_tsv("./codes-phase.inp")

phase_df <- 
    tibble(subject = map_chr(best_list, 1),
	   locus = paste(loci_df$locus, collapse = " "),
	   h1 = trimws(map_chr(best_list, 2)),
	   h2 = trimws(map_chr(best_list, 3))) %>%
    mutate(subject = sub("0 #", "", subject)) %>%
    separate_rows(locus, h1, h2, sep = " ") %>%
    mutate(uncertain = as.integer(grepl("\\(", h1))) %>%
    mutate_at(vars(h1, h2), . %>% gsub("\\(|\\)", "", .)) %>%
    gather(hap, allele, h1:h2) %>%
    mutate(hap = gsub("h", "", hap),
	   allele = as.integer(allele)) %>%
    select(subject, locus, hap, allele, uncertain) %>%
    left_join(loci_df, by = "locus") %>%
    select(subject, locus, pos, hap, allele, uncertain) %>%
    arrange(subject, pos, hap) %>%
    left_join(hla_codes, by = c("subject", "locus", "allele" = "code")) %>%
    mutate(allele = ifelse(!is.na(allele.y), allele.y, allele)) %>%
    select(subject, locus, pos, hap, allele, uncertain) %>%
    distinct()

hla_qtls <- 
    "../eqtl_mapping/transcriptomemapping/hla_personalized/2-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_best == 1L) %>%
    select(phen_id, rank, var_id) %>%
    inner_join(select(gencode_hla, gene_name, gene_id), by = c("phen_id" = "gene_id")) %>%
    select(locus = gene_name, rank, rsid = var_id)

hla_df <- phase_df %>%
    filter(grepl("HLA", locus))

qtl_df <- phase_df %>%
    filter(!grepl("HLA", locus)) %>%
    rename(rsid = locus) %>%
    left_join(hla_qtls, by = "rsid")

gene_snp_df <- left_join(hla_df, qtl_df, by = c("subject", "locus", "hap"),
			 suffix = c("_gene", "_snp"))

write_tsv(gene_snp_df, "./phase_hla_haps_snps.tsv")
