#!/bin/bash

JOBS=12
parallel=/home/vitor/parallel
samples=../../data/sample_info/samples_phase3_ena.txt
write_custom_index=../../../imgt_index/write_genotyped_alleles.R

$parallel --gnu -j $JOBS ./map_and_quantify_round1.sh {} :::: $samples

Rscript process_quants.R 1

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./sample_indices

$parallel --gnu -j $JOBS ./map_and_quantify_round2.sh {} :::: $samples

Rscript process_quants.R 2

Rscript make_expression_matrix.R
