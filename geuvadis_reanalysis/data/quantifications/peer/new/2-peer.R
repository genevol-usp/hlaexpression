library(magrittr)
library(peer)

k <- 10L

phenotypes <- 
  data.table::fread("./geuvadis_fpkms.csv") %>%
  as.data.frame() %>%
  `rownames<-`(.$subject) %>%
  .[-1] %>%
  data.matrix()

model <- PEER()
PEER_setPhenoMean(model, phenotypes)
PEER_setNk(model, k)
PEER_setNmax_iterations(model, 1e6)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)

peer_residuals <- PEER_getResiduals(model)
colnames(peer_residuals) <- colnames(phenotypes)

residuals_df <- tibble::as_tibble(peer_residuals) %>%
  tibble::add_column(subject = rownames(phenotypes), .before = 1)

data.table::fwrite(residuals_df, paste0("./residuals_", k, ".tsv"),
		   sep = "\t", quote = FALSE)
