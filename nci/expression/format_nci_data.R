library(tidyverse)

nci <- 
  "../data/HLA-A, -B, -C expression levels for Pat.xlsx" %>%
  readxl::read_excel() %>%
  select(3, 5:14)

mRNA <-
  select(nci, subject = 1, ends_with("mRNA")) %>%
  gather(locus, mRNA, -1) %>%
  mutate(locus = sub(" mRNA", "", locus),
	 mRNA = round(as.numeric(mRNA), 3)) %>%
  arrange(subject, locus)

c_surface <- select(nci, subject = 1, c_surface = `HLA-C Surface`) 

genos <-
  select(nci, subject = 1, starts_with("Class")) %>%
  gather(variable, allele, -1) %>%
  extract(variable, c("locus", "hap"), "Class 1([ABC]) [ABC](\\d)") %>%
  arrange(subject, locus, hap) %>%
  mutate(allele = ifelse(nchar(allele) == 5, 
			 sub("^(\\d{4})(\\d)$", "\\10\\2", allele), allele)) %>%
  separate(allele, c("f1", "f2", "f3"), c(2, 4)) %>% 
  unite(allele, f1, f2, f3, sep = ":") %>% 
  mutate(allele = gsub(":$|:00:$", "", allele),
	 allele = paste0(locus, "*", allele),
	 locus = paste0("HLA-", locus))

genos$allele[genos$subject == '66K00241' & genos$allele == 'A*03:01:01'] <- "A*30:01:01"

nci_merge <- 
  left_join(genos, mRNA, by = c("subject", "locus")) %>%
  left_join(c_surface, by = c("subject"))

write_tsv(nci_merge, "../data/nci_expression.tsv")
