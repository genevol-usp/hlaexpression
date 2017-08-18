#!/bin/bash

JOBS=10
parallel=/home/vitor/parallel
samples=$(echo sample_{01..50})

mismatch=0.1

$parallel --gnu -j $JOBS ./map_and_quantify_CHR.sh {} $mismatch ::: $samples

Rscript process_quants.R CHR_$mismatch
