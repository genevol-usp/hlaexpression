devtools::load_all("~/hlaseqlib")
library(tidyverse)

locus <- commandArgs(TRUE)[1]
gen_file <- paste0("~/IMGTHLA/alignments/", locus, "_gen.txt")

if (grepl("DRB\\d+", locus)) { 
    locus_nuc <- "DRB"
} else {
    locus_nuc <- locus
}

nuc_file <- paste0("~/IMGTHLA/alignments/", locus_nuc, "_nuc.txt")

hla_df <- hla_read_alignment(nuc_file) %>%
   mutate(gene = sub("^([^*]+).+$", "\\1", allele)) %>%
   filter(gene == locus) %>%
   select(-gene)

if (all(grepl("\\*", hla_df$cds))) {
    stop(paste("no complete sequence for locus", locus))
}

if (nrow(hla_df) == 1L || all(!grepl("\\*", hla_df$cds))) {
    
    final_df <- hla_df

} else {

    distmatrix <- make_dist_matrix(hla_df)

    closest_allele_df <- make_closest_allele_df(distmatrix)

    closest_allele_df$id <- closest_allele_df %>% group_indices(inc_allele)

    closest_allele_df_step2 <-
	bind_rows(select(closest_allele_df, id, allele = inc_allele),
		  select(closest_allele_df, id, allele = closest)) %>%
	distinct() %>%
	left_join(hla_df, by = "allele") %>%
	split(.$id) %>%
	map(make_dist_matrix) %>%
	map(make_closest_allele_df) %>%
	bind_rows()

    closest_within_type <- closest_allele_df_step2 %>%
	find_closest_within_type()

    inferred_df <- closest_within_type %>%
	left_join(hla_df, by = c("inc_allele" = "allele")) %>%
	left_join(hla_df, by = c("closest" = "allele")) %>%
	mutate(cds = purrr::map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
	select(allele = inc_allele, cds)

    final_df <- hla_df %>%
	filter(!grepl("\\*", cds)) %>%
	bind_rows(inferred_df) %>%
	arrange(allele)
}

if (!file.exists(gen_file)) {

    out_df <- final_df %>%
	mutate(cds = hla_format_sequence(cds)) %>%
	rename(transcript = cds) %>%
	mutate(allele3f = hla_trimnames(allele, 3)) %>%
	distinct(allele3f, transcript, .keep_all = TRUE) %>%
	group_by(allele3f) %>%
	mutate(n = n()) %>%
	ungroup() %>%
	mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
	select(allele, transcript) %>%
	arrange(allele)

} else {

    utr_df <- hla_read_utr(gen_file) %>%
	full_join(final_df, by = "allele") %>%
	mutate_at(vars(utr5, utr3), 
		  ~ifelse(is.na(.), paste(rep("*", 100L), collapse = ""), .)) 
	
    hla_df_utr <- utr_df %>%
	unite(cds, utr5, cds, utr3, sep = "")

    utr_df <- utr_df %>%
	unite(utr, utr5, utr3, sep = "|") %>%
	select(allele, utr)

    if (all(!grepl("\\*", utr_df$utr))) {
	
	final_utr_df <- utr_df %>%
	    separate(utr, c("utr5", "utr3"), sep = "\\|")
    
    } else {

	distmatrix_utr <- make_dist_matrix(hla_df_utr)

	closest_allele_utr <- make_closest_allele_df(distmatrix_utr)

	closest_within_type_utr <- closest_allele_utr %>%
	    find_closest_within_type()

	inferred_utr_df <- closest_within_type_utr %>%
	    left_join(utr_df, by = c("inc_allele" = "allele")) %>% 
	    left_join(utr_df, by = c("closest" = "allele")) %>%
	    mutate(utr = map2_chr(utr.x, utr.y, hla_attribute_seq)) %>%
	    select(allele = inc_allele, utr)

	final_utr_df <- utr_df %>%
	    filter(!grepl("\\*", utr)) %>%
	    bind_rows(inferred_utr_df) %>%
	    arrange(allele) %>%
	    separate(utr, c("utr5", "utr3"), sep = "\\|")
    }

    out_df <- 
	left_join(final_df, final_utr_df, by = "allele") %>%
	select(allele, utr5, cds, utr3) %>%
	mutate_at(vars(utr5, cds, utr3), hla_format_sequence) %>%
	mutate(utr5 = substring(utr5, nchar(utr5)-74L, nchar(utr5)),
	       utr3 = substring(utr3, 1, 75)) %>%
	unite(transcript, utr5, cds, utr3, sep = "") %>%
	mutate(allele3f = hla_trimnames(allele, 3)) %>%
	distinct(allele3f, transcript, .keep_all = TRUE) %>%
	group_by(allele3f) %>%
	mutate(n = n()) %>%
	ungroup() %>%
	mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
	select(allele, transcript) %>%
	arrange(allele)
}

write_tsv(out_df, paste0("./index_", locus, ".tsv"))
