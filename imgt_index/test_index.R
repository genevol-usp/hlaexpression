devtools::load_all("~/hlaseqlib")
    
make_closest_allele_df <- function(distmatrix, inc_alleles, comp_alleles) {

    distmatrix[inc_alleles, comp_alleles] %>%
    apply(1, function(x) which(x == min(x, na.rm = TRUE))) %>%
    purrr::map(~tibble::tibble(closest = names(.))) %>%
    dplyr::bind_rows(.id = "inc_allele") %>% 
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

locus <- "A"
gen_file <- paste0("~/IMGTHLA/alignments/", locus, "_gen.txt")
nuc_file <- paste0("~/IMGTHLA/alignments/", locus, "_nuc.txt")

hla_df <- hla_read_alignment("~/IMGTHLA/alignments/A_nuc.txt") 

hla_df$cds <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
    apply(1, . %>% paste(collapse = ""))

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

cds_common <- hla_df_cds_common$cds
names(cds_common) <- hla_df_cds_common$allele

distmatrix <- stringdist::stringdistmatrix(cds_common, cds_common,
					   method = "hamming", 
					   useNames = "names")

diag(distmatrix) <- NA

min_dist_df <-
    make_closest_allele_df(distmatrix, incomplete_alleles, complete_alleles)

inferred_df <- min_dist_df %>%
    dplyr::left_join(hla_df, by = c("inc_allele" = "allele")) %>%
    dplyr::left_join(hla_df, by = c("closest" = "allele")) %>%
    dplyr::mutate(cds = purrr::map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
    dplyr::select(allele = inc_allele, cds)

final_df <- hla_df %>%
    dplyr::filter(allele %in% complete_alleles) %>%
    dplyr::bind_rows(inferred_df) %>%
    dplyr::arrange(allele)

utr_df <- hla_read_utr(gen_file) %>%
    dplyr::full_join(tibble::tibble(allele = final_df$allele), by = "allele") %>%
    dplyr::mutate_at(dplyr::vars(utr5, utr3), 
		     ~ifelse(is.na(.), paste(rep("*", 100L), collapse = ""), .)) %>%
    tidyr::unite(utr, utr5, utr3, sep = "|")

complete_utr_alleles <- utr_df$allele[!grepl("\\*", utr_df$utr)]
incomplete_utr_alleles <- utr_df$allele[grepl("\\*", utr_df$utr)]

min_dist_utr_df <- 
    make_closest_allele_df(distmatrix, incomplete_utr_alleles, complete_utr_alleles) 

inferred_utr_df <- min_dist_utr_df %>%
    dplyr::left_join(utr_df, by = c("inc_allele" = "allele")) %>%
    dplyr::left_join(utr_df, by = c("closest" = "allele")) %>%
    dplyr::mutate(utr = purrr::map2_chr(utr.x, utr.y, hla_attribute_seq)) %>%
    dplyr::select(allele = inc_allele, utr)

final_utr_df <- utr_df %>%
    dplyr::filter(allele %in% complete_utr_alleles) %>%
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

