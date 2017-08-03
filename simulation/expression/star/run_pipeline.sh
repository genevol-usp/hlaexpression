#!/bin/bash

JOBS=10
parallel=/home/vitor/parallel
samples=$(echo sample_{01..50})

write_custom_index=../../../geuvadis_reanalysis/expression/kallisto/index/write_genotyped_alleles.R

$parallel --gnu -j $JOBS ./map_and_quantify_round1.sh {} ::: $samples

Rscript process_quants.R 1

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./sample_indices

$parallel --gnu -j $JOBS ./map_and_quantify_round2.sh {} ::: $samples

Rscript process_quants.R 2
