#!/bin/bash

JOBS=50
parallel=/home/vitor/parallel

$parallel --gnu -j $JOBS ./kallisto_quant_CHR.sh {} ::: $(echo sample_{01..50})

Rscript process_kallisto.R CHR

$parallel --gnu -j $JOBS ./kallisto_quant_ALL.sh {} ::: $(echo sample_{01..50})

Rscript process_kallisto.R ALL
