#!/bin/bash

JOBS=12
parallel=/home/vitor/parallel
samples=$(echo sample_{01..50})

$parallel --gnu -j $JOBS ./map_and_quantify_CHR.sh {} ::: $samples

Rscript process_quants.R CHR
