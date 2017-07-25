#!/bin/bash

JOBS=30
parallel=/home/vitor/parallel
samples=../../data/sample_info/samples_phase3_ena.txt

$parallel --gnu -j $JOBS ./kallisto_quant_CHR.sh {} :::: $samples

Rscript process_kallisto.R CHR

$parallel --gnu -j $JOBS ./kallisto_quant_ALL.sh {} :::: $samples

Rscript process_kallisto.R ALL
