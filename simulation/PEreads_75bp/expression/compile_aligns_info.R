library(tidyverse)

samples <- sprintf("sample_%02d", 1:50)

not_aligned_star_pri <-
    file.path("./star/mappings_PRI", samples, "reads_not_aligned_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

not_aligned_star_imgt <-
    file.path("./star/mappings_2", samples, "reads_not_aligned_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

not_aligned_kallisto_pri <-
    file.path("./kallisto/quantifications_PRI", samples, "reads_not_aligned_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

not_aligned_kallisto_imgt <-
    file.path("./kallisto/quantifications_2", samples, "reads_not_aligned_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

list(star_pri = not_aligned_star_pri, 
     star_imgt = not_aligned_star_imgt,
     kallisto_pri = not_aligned_kallisto_pri,
     kallisto_imgt = not_aligned_kallisto_imgt) %>%
    bind_rows(.id = "aligner") %>%
    group_by(gene_read, aligner) %>%
    summarize(na = mean(na)) %>%
    spread(aligner, na) %>%
    write_tsv("./reads_not_aligned_hla.tsv")

falseneg_star_pri <-
    file.path("./star/mappings_PRI", samples, "reads_false_neg_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falseneg_star_imgt <- 
    file.path("./star/mappings_2", samples, "reads_false_neg_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falseneg_kallisto_pri <-
    file.path("./kallisto/quantifications_PRI", samples, "reads_false_neg_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falseneg_kallisto_imgt <-
    file.path("./kallisto/quantifications_2", samples, "reads_false_neg_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

list(star_pri = falseneg_star_pri,
     star_imgt = falseneg_star_imgt,
     kallisto_pri = falseneg_kallisto_pri,
     kallisto_imgt = falseneg_kallisto_imgt) %>%
    bind_rows(.id = "aligner") %>%
    group_by(gene_read, aligner) %>%
    summarize(perc = mean(perc)) %>%
    spread(aligner, perc) %>%
    write_tsv("./reads_false_neg_hla.tsv")

falsepos_star_pri <-
    file.path("./star/mappings_PRI", samples, "reads_false_pos_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falsepos_star_imgt <-
    file.path("./star/mappings_2", samples, "reads_false_pos_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falsepos_kallisto_pri <-
    file.path("./kallisto/quantifications_PRI", samples, "reads_false_pos_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

falsepos_kallisto_imgt <-
    file.path("./kallisto/quantifications_2", samples, "reads_false_pos_hla.tsv") %>%
    setNames(samples) %>%
    plyr::ldply(read_tsv, .id = "subject")

list(star_pri = falsepos_star_pri,
     star_imgt = falsepos_star_imgt,
     kallisto_pri = falsepos_kallisto_pri,
     kallisto_imgt = falsepos_kallisto_imgt) %>%
    bind_rows(.id = "aligner") %>%
    group_by(gene_ref, aligner) %>%
    summarize(perc = mean(perc)) %>%
    spread(aligner, perc) %>%
    write_tsv("./reads_false_pos_hla.tsv")
