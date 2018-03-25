library(tidyverse)

D <- read_delim("./enrichment_qtl_in_tf.txt", col_names = FALSE, delim = " ")

fisher.test(matrix(c(D$X1, D$X2, round(D$X3), D$X2), ncol=2))
