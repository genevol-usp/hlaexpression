#!/bin/bash

parallel=/home/vitor/parallel

JOBS=30

$parallel --gnu -j $JOBS ./kallisto_quant_CHR.sh {} :::: samples_phase3_ena.txt

#Rscript process_kallisto.R CHR

#$parallel --gnu -j $JOBS ./kallisto_quant_ALL.sh {} :::: samples_phase3_ena.txt

#Rscript process_kallisto.R ALL
