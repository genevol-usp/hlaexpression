library(data.table)
library(Biostrings)

imgt <- readDNAStringSet("../../geuvadis_reanalysis/expression/kallisto/index/imgt_index.fa")
imgt <- imgt[grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)\\*", names(imgt))]

bed <- data.table(chr = names(imgt),
		  start = 0L,
		  end = width(imgt))

fwrite(bed, "./imgt.bed", sep = "\t", quote = FALSE, col.names = FALSE)
