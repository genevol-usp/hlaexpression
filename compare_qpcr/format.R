library(tidyverse)
library(readxl)

a_alleles <- 
    tribble(~g, ~lineage,
	    1, "A*03",
	    2, "A*32",
	    3, "A*33",
	    4, "A*31",
	    5, "A*25",
	    6, "A*30",
	    7, "A*02", 
	    8, "A*29",
	    9, "A*11",
	    10, "A*26",
	    11, "A*68",
	    12, "A*01",
	    13, "A*23",
	    14, "A*24")

b_alleles <- 
    tribble(~g, ~lineage,
	    1, "B*52",
	    2, "B*15",
	    3, "B*51",
	    4, "B*18",
	    5, "B*08",
	    6, "B*57",
	    7, "B*07",
	    8, "B*40",
	    9, "B*14",
	    10, "B*44",
	    11, "B*55",
	    12, "B*27",
	    13, "B*49",
	    14, "B*38",
	    15, "B*35",
	    16, "B*39",
	    17, "B*37",
	    18, "B*13")



hla_a <- read_csv("./HLA_A.csv", col_names = FALSE) %>% 
    arrange(X1) %>% 
    mutate(X3 = c(0.5, X1[1:(length(X1) - 1)]),
	   g = cumsum(X1 - X3 >= 0.5) + 1) %>%
    left_join(a_alleles, by = "g") %>%
    select(lineage, mrna = X2)

hla_b <- read_csv("./HLA_B.csv", col_names = FALSE) %>%
    arrange(X1) %>%
    mutate(X3 = c(0.5, X1[1:(length(X1) - 1)]),
	   g = cumsum(X1 - X3 >= 0.4) + 1) %>%
    left_join(b_alleles, by = "g") %>%
    select(lineage, mrna = X2)

hla_c <- read_excel("./HLA_C.xlsx") %>%
    select(subject = 1, allele1 = 3, allele2 = 4, mrna = 2) %>%
    mutate(mrna = as.numeric(mrna)) %>%
    filter(!is.na(mrna)) %>%
    gather(i, allele, allele1:allele2) %>%
    mutate(lineage = sub("^(\\d+).*$", "\\1", allele),
           lineage = paste0("C*", lineage)) %>%
    select(subject, lineage, allele, mrna) %>%
    arrange(subject, lineage, allele) 

write_tsv(hla_a, "HLA_A_formatted.tsv")
write_tsv(hla_b, "HLA_B_formatted.tsv")
write_tsv(hla_c, "HLA_C_formatted.tsv")
