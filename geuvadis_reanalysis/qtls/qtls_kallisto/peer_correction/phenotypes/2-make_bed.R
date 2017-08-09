devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

write_bed <- function(peer_out) { 
 
  peerk <- sub("^residuals_(\\d+)\\.tsv$", "\\1", basename(peer_out))
  
  bed <-
    data.table::fread(peer_out) %>%
    gather(gene_id, resid, -subject) %>%
    spread(subject, resid) %>% 
    inner_join(gencode_chr_gene, by = "gene_id") %>% 
    filter(chr %in% 1:22) %>% 
    mutate(chr = as.integer(chr), id = gene_id) %>%
    select(`#chr` = chr, start, end, id, gid = gene_id, strd = strand, 
	   starts_with("HG"), starts_with("NA")) %>%
    arrange(`#chr`, start, end)

  bed_out <- sprintf("./phenotypes_eur_%s.bed", peerk)
  bed_out_gz <- paste0(bed_out, ".gz")
  write_tsv(bed, bed_out)
  system(paste("bgzip", bed_out, "&& tabix -p bed", bed_out_gz))
}

peer_files <- list.files(".", pattern = "residuals_\\d+\\.tsv$", full.names = TRUE)

parallel::mclapply(peer_files, write_bed, mc.cores = length(peer_files))
