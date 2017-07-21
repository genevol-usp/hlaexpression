#!/bin/bash

JOBS=60
parallel=/home/vitor/parallel

#$parallel --gnu -j $JOBS ./kallisto_quant_round1.sh {} :::: samples_phase3_ena.txt

#Rscript process_kallisto.R 1

Rscript ./index/write_genotyped_alleles.R\
  ./quantifications_1/processed_quant.tsv ./index/sample_indices

$parallel --gnu -j $JOBS ./kallisto_quant_round2.sh {} :::: samples_phase3_ena.txt

Rscript process_kallisto.R 2

Rscript make_kallisto_matrix.R
