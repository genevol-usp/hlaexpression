devtools::load_all("/home/vitor/genomicRutils")
library(tidyverse)

write_bed <- function(peer_out) { 
 
  peerk <- sub("^.*residuals_(\\d+)\\.rds$", "\\1", peer_out)
  peer_mat <- readRDS(peer_out) 
  subjects <- rownames(peer_mat)
  
  bed <-
    as_tibble(peer_mat) %>%
    add_column(subject = subjects, .before = 1) %>%
    gather(id, resid, -subject) %>% 
    spread(subject, resid) %>% 
    inner_join(gencode_chr_gene, by = c("id" = "gene_id")) %>% 
    filter(chr %in% 1:22) %>% 
    mutate(chr = as.integer(chr), gid = id) %>%
    select(`#chr` = chr, start, end, id, gid, strd = strand, starts_with("HG"), starts_with("NA")) %>%
    arrange(`#chr`, start, end)

  bed_out <- sprintf("./phenotypes/phenotypes_eur_%s.bed", peerk)
  bed_out_gz <- paste0(bed_out, ".gz")
  write_tsv(bed, bed_out)
  system(paste("bgzip", bed_out, "&& tabix -p bed", bed_out_gz))
}

peer_files <- list.files("../../peer", pattern = "\\.rds$", full.names = TRUE)

parallel::mclapply(peer_files, write_bed, mc.cores = length(peer_files))
