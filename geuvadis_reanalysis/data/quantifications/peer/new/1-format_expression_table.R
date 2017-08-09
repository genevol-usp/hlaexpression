devtools::load_all("/home/vitor/hlaseqlib")
library(data.table)

autosomes <- 
  get_gencode_coords("/home/vitor/gencode_data/gencode.v19.annotation.gtf.gz", 
		     feature = "gene") %>%
  dplyr::filter(chr %in% 1:22) %>%
  dplyr::pull(gene_id)

quant <- 
  fread("../../geuvadis_667samples_gencord19_rpkm_358samples.txt", drop = 2)
setnames(quant, 1, "gene_id")

quant <- quant[gene_id %in% autosomes]

quant_l <- melt(quant, id = 1, measure.vars = 2:ncol(quant), 
	      variable.name = "subject", value.name = "fpkm") 

expressed_genes <- quant_l[, .(mean(fpkm > 0)), by = .(gene_id)][V1 >= 0.9]

quant_l <- quant_l[gene_id %in% expressed_genes$gene_id]

quant_w <- dcast(quant_l, subject ~ gene_id, value.var = "fpkm") 

quant_w[, subject := sub("^([^\\.]+).+$", "\\1", subject)]

fwrite(quant_w, "./geuvadis_fpkms.csv")
