#!/bin/bash

JOBS=50
parallel=/home/vitor/parallel

write_custom_idx=../../../geuvadis_reanalysis/expression/kallisto/index/write_genotyped_alleles.R

$parallel --gnu -j $JOBS ./kallisto_quant_round1.sh {} ::: $(echo sample_{01..50})

Rscript process_kallisto.R 1

Rscript $write_custom_idx ./quantifications_1/processed_quant.tsv ./index/sample_indices

$parallel --gnu -j $JOBS ./kallisto_quant_round2.sh {} ::: $(echo sample_{01..50})

Rscript process_kallisto.R 2
