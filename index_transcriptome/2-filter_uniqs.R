library(Biostrings)

index <- readDNAStringSet("./gencode.v25.PRI.transcripts.fa")

writeXStringSet(unique(index), "./gencode.v25.PRI.uniqTranscripts.fa")
