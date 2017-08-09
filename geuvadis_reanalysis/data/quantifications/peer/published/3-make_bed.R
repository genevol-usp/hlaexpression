devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gencode12 <- 
  "/home/vitor/gencode_data/gencode.v12.annotation.gtf.gz" %>%
  get_gencode_coords(feature = "gene")

bed <-
  data.table::fread("./residuals_10.tsv") %>%
  gather(gene_id, resid, -subject) %>%
  group_by(gene_id) %>%
  mutate(resid = GenABEL::rntransform(resid)) %>%
  ungroup() %>%
  spread(subject, resid) %>% 
  inner_join(gencode12, by = "gene_id") %>% 
  filter(chr %in% 1:22) %>% 
  mutate(chr = as.integer(chr), id = gene_id) %>%
  select(`#chr` = chr, start, end, id, gid = gene_id, strd = strand, 
	 starts_with("HG"), starts_with("NA")) %>%
  arrange(`#chr`, start, end)

bed_out <- "./phenotypes_eur_10.bed"
bed_out_gz <- paste0(bed_out, ".gz")
write_tsv(bed, bed_out)
system(paste("bgzip", bed_out, "&& tabix -p bed", bed_out_gz))
