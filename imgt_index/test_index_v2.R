devtools::load_all("~/hlaseqlib")

# Functions
make_dist_matrix <- function(hla_df) {

    complete_alleles <- hla_df$allele[!grepl("\\*", hla_df$cds)]
    incomplete_alleles <- hla_df$allele[grepl("\\*", hla_df$cds)]

    cds_sequenced <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
	apply(2, function(x) !any(x == "*"))

    run <- rle(cds_sequenced)

    ends <- cumsum(run$lengths)
    starts <- ends - run$lengths + 1L

    run_df <- 
	tibble::tibble(value = run$values, 
		       start = starts, 
		       end = ends) %>%
	dplyr::filter(value == TRUE) %>%
	dplyr::mutate(l = purrr::map2_int(start, end, ~length(.x:.y))) %>%
	dplyr::filter(l == max(l))

    hla_df_cds_common <- hla_df %>%
	dplyr::mutate(cds = substring(cds, run_df$start, run_df$end))

    hla_df_cds_common_complete <- hla_df_cds_common %>%
	dplyr::filter(allele %in% complete_alleles)

    hla_df_cds_common_incomplete <- hla_df_cds_common %>%
	dplyr::filter(allele %in% incomplete_alleles)

    cds_common_complete <- hla_df_cds_common_complete$cds
    names(cds_common_complete) <- hla_df_cds_common_complete$allele
    
    cds_common_incomplete <- hla_df_cds_common_incomplete$cds
    names(cds_common_incomplete) <- hla_df_cds_common_incomplete$allele

    stringdist::stringdistmatrix(cds_common_incomplete, cds_common_complete,
				 method = "hamming", useNames = "names")
}
    
make_closest_allele_df <- function(distmatrix) {

    cnames <- colnames(distmatrix)
    
    distmatrix %>%
	split(seq_len(nrow(distmatrix))) %>%
	setNames(rownames(distmatrix)) %>%
	lapply(function(x) which(x == min(x, na.rm = TRUE))) %>% 
	purrr::map(~tibble::tibble(closest = cnames[.])) %>%
	dplyr::bind_rows(.id = "inc_allele")
}

find_closest_within_type <- function(closest_df) {

    closest_df %>% 
    dplyr::mutate(`1` = hla_trimnames(inc_allele, 1) == hla_trimnames(closest, 1),
		  `2` = hla_trimnames(inc_allele, 2) == hla_trimnames(closest, 2),
		  `3` = hla_trimnames(inc_allele, 3) == hla_trimnames(closest, 3),
		  `4` = hla_trimnames(inc_allele, 4) == hla_trimnames(closest, 4)) %>%
    tidyr::gather(field, value, `1`:`4`) %>%
    dplyr::group_by(inc_allele) %>%
    dplyr::filter(all(value == FALSE) | value == TRUE) %>%
    dplyr::slice(which.max(field)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(inc_allele) %>%
    dplyr::select(inc_allele, closest)
}


# Analysis
locus <- "A"
gen_file <- paste0("~/IMGTHLA/alignments/", locus, "_gen.txt")
nuc_file <- paste0("~/IMGTHLA/alignments/", locus, "_nuc.txt")

hla_df <- hla_read_alignment("~/IMGTHLA/alignments/A_nuc.txt") 

distmatrix <- make_dist_matrix(hla_df)

closest_allele_df <- make_closest_allele_df(distmatrix)

closest_allele_df$id <- closest_allele_df %>% dplyr::group_indices(inc_allele)

closest_allele_df_step2 <-
    dplyr::bind_rows(dplyr::select(closest_allele_df, id, allele = inc_allele),
		     dplyr::select(closest_allele_df, id, allele = closest)) %>%
    dplyr::distinct() %>%
    dplyr::left_join(hla_df, by = "allele") %>%
    split(.$id) %>%
    purrr::map(make_dist_matrix) %>%
    purrr::map(make_closest_allele_df) %>%
    dplyr::bind_rows()

closest_within_type <- closest_allele_df_step2 %>%
    find_closest_within_type()

inferred_df <- closest_within_type %>%
    dplyr::left_join(hla_df, by = c("inc_allele" = "allele")) %>%
    dplyr::left_join(hla_df, by = c("closest" = "allele")) %>%
    dplyr::mutate(cds = purrr::map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
    dplyr::select(allele = inc_allele, cds)

final_df <- hla_df %>%
    dplyr::filter(!grepl("\\*", cds)) %>%
    dplyr::bind_rows(inferred_df) %>%
    dplyr::arrange(allele)

utr_df <- hla_read_utr(gen_file) %>%
    dplyr::full_join(final_df, by = "allele") %>%
    dplyr::mutate_at(dplyr::vars(utr5, utr3), 
		     ~ifelse(is.na(.), paste(rep("*", 100L), collapse = ""), .)) 
    
hla_df_utr <- utr_df %>%
    tidyr::unite(cds, utr5, cds, utr3, sep = "")

utr_df <- utr_df %>%
    tidyr::unite(utr, utr5, utr3, sep = "|") %>%
    dplyr::select(allele, utr)

distmatrix_utr <- make_dist_matrix(hla_df_utr)

closest_allele_utr <- make_closest_allele_df(distmatrix_utr)

closest_within_type_utr <- closest_allele_utr %>%
    find_closest_within_type()

inferred_utr_df <- closest_within_type_utr %>%
    dplyr::left_join(utr_df, by = c("inc_allele" = "allele")) %>% 
    dplyr::left_join(utr_df, by = c("closest" = "allele")) %>%
    dplyr::mutate(utr = purrr::map2_chr(utr.x, utr.y, hla_attribute_seq)) %>%
    dplyr::select(allele = inc_allele, utr)

final_utr_df <- utr_df %>%
    dplyr::filter(!grepl("\\*", utr)) %>%
    dplyr::bind_rows(inferred_utr_df) %>%
    dplyr::arrange(allele) %>%
    tidyr::separate(utr, c("utr5", "utr3"), sep = "\\|")

out_df <- 
    dplyr::left_join(final_df, final_utr_df, by = "allele") %>%
    dplyr::select(allele, utr5, cds, utr3) %>%
    dplyr::mutate_at(dplyr::vars(utr5, cds, utr3), hla_format_sequence) %>%
    dplyr::mutate(utr5 = substring(utr5, nchar(utr5)-74L, nchar(utr5)),
		  utr3 = substring(utr3, 1, 75)) %>%
    tidyr::unite(transcript, utr5, cds, utr3, sep = "") %>%
    dplyr::mutate(allele3f = hla_trimnames(allele, 3)) %>%
    dplyr::distinct(allele3f, transcript, .keep_all = TRUE) %>%
    dplyr::group_by(allele3f) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
    dplyr::select(allele, transcript) %>%
    dplyr::arrange(allele)

