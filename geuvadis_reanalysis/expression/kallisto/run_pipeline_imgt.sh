#!/bin/bash

JOBS=60
parallel=/home/vitor/parallel
samples=../../data/sample_info/samples_phase3_ena.txt
write_custom_index=./index/write_genotyped_alleles.R 

$parallel --gnu -j $JOBS ./kallisto_quant_round1.sh {} :::: $samples

Rscript process_kallisto.R 1

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./index/sample_indices

$parallel --gnu -j $JOBS ./kallisto_quant_round2.sh {} :::: $samples

Rscript process_kallisto.R 2

Rscript make_kallisto_matrix.R
