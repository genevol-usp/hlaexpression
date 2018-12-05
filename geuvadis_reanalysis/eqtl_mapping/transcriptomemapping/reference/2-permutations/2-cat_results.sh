#!/bin/bash

runFDR=/home/vitor/QTLtools/script/runFDR_cis.R
PC=60
OUT=./results/permutations_$PC

cat ${OUT}_*.txt | gzip -c > $OUT.txt.gz
rm ${OUT}_*.txt

Rscript $runFDR $OUT.txt.gz 0.05 $OUT
