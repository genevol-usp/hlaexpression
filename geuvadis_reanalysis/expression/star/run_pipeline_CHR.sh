#!/bin/bash

JOBS=12
parallel=/home/vitor/parallel
samples=../../data/sample_info/samples_phase3_ena.txt

$parallel --gnu -j $JOBS ./map_and_quantify_CHR.sh {} :::: $samples

Rscript process_quants.R CHR
