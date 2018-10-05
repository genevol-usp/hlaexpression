library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib")

#info <- "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB3366&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"
#download.file(info, destfile = "./run_info.txt")

fastqs <- read_tsv("./run_info.txt") %>% 
    extract(submitted_ftp, "run_id", "fastq/([^;]+)") %>%
    extract(run_id, "subject_id", "([^.]+)", remove = FALSE) %>%
    mutate(run_id = sub("_\\d\\.fastq\\.gz$", "", run_id)) %>%
    select(subject_id, run_id, ena_id = run_accession, fastq_ftp)
 
final_df <- geuvadis_info %>%
    select(subject_id = name, run_id = assay_name, ena_id)

all_inds <- 
    read_tsv("../../data/quantifications/GD660.GeneQuantRPKM.txt.gz", n_max = 0) %>%
    select(-(1:4)) %>%
    do(tibble(run_id = names(.))) %>%
    mutate(subject_id = sub("^([^.]+).+$", "\\1", run_id)) %>%
    left_join(fastqs, by = c("subject_id", "run_id")) %>%
    mutate(final_set = as.integer(run_id %in% final_df$run_id)) %>%
    arrange(ena_id)

write_tsv(select(all_inds, subject_id, ena_id), "./sample_info.tsv")

replicates <- all_inds %>% 
    group_by(subject_id) %>% 
    filter(n() == 2) %>% 
    ungroup() %>%
    arrange(subject_id)

ftp_df <- replicates %>%
    filter(final_set == 0) %>%
    left_join(select(geuvadis_info, subject_id = name, kgp_phase3, pop),
	      by = "subject_id") %>%
    filter(pop != "YRI", kgp_phase3 == 1) 

ftp_urls <- ftp_df %>%
    separate_rows(fastq_ftp, sep = ";") %>%
    pull(fastq_ftp)

writeLines(ftp_urls, "./ftp_urls.txt")

ftp_df %>% 
    arrange(ena_id) %>%
    pull(ena_id) %>%
    writeLines("./sample_ids.tsv")
