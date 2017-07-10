#!/bin/bash

JOBS=36
parallel=/home/vitor/parallel

$parallel --gnu -j $JOBS ./kallisto_quant_round1.sh {} :::: nci_samples.txt

Rscript process_kallisto.R 1

Rscript ../index/write_genotyped_alleles.R\
  ./quantifications_1/processed_quant.tsv ./sample_indices

$parallel --gnu -j $JOBS ./kallisto_quant_round2.sh {} :::: ./nci_samples.txt

Rscript process_kallisto.R 2
