devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <-   
  geuvadis_info %>%
  filter(kgp_phase3 == 1L & pop != "YRI") %>%
  pull(assay_name)

GD462 <- 
  data.table::fread("zcat < ../GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")

GD660 <- data.table::fread("zcat < ../GD660.GeneQuantRPKM.txt.gz") %>%
  filter(Gene_Symbol %in% GD462$Gene_Symbol) %>%
  select(gene_id = Gene_Symbol, samples) %>%
  setNames(sub("^([^\\.]+).*$", "\\1", names(.))) %>%
  gather(subject, value, -gene_id) %>%
  spread(gene_id, value)

data.table::fwrite(GD660, "./geuvadis_rpkms.tsv", sep = "\t", quote = FALSE)
