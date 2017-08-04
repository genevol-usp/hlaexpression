library(magrittr)
library(peer)

k <- 10 

phenotypes <- 
  data.table::fread("../../../../expression/kallisto/kallisto_gene_expressed90%_fpkm.csv") %>%
  as.data.frame() %>%
  `rownames<-`(.$subject) %>%
  dplyr::select(-subject) %>%
  data.matrix()

covs <- 
  data.table::fread("./covs.csv") %>%
  as.data.frame() %>%
  `rownames<-`(.$subject) %>%
  dplyr::select(-subject) %>%
  data.matrix() %>%
  .[rownames(phenotypes), ]
  
model <- PEER()
PEER_setPhenoMean(model, phenotypes)
PEER_setCovariates(model, covs)
PEER_setNk(model, k)
PEER_setNmax_iterations(model, 1e6)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)

peer_residuals <- PEER_getResiduals(model)
dimnames(peer_residuals) <- dimnames(phenotypes)

resid_mean_added <- 
  apply(phenotypes, 2, mean) %>%
  sweep(peer_residuals, 2, ., `+`)

saveRDS(resid_mean_added, paste0("./residuals/residuals_", k, ".rds"))
