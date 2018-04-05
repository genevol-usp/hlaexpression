#!/bin/bash

QTLtools_dir=/home/vitor/QTLtools
QTLtools=$QTLtools_dir/QTLtools_1.1_Ubuntu16.04_x86_64

for PC in $(seq 0 5 20; seq 30 10 100)
do
  OUT=./results/permutations_$PC
  cat ${OUT}_*.txt | gzip -c > $OUT.txt.gz
  rm ${OUT}_*.txt

  Rscript $QTLtools_dir/script/runFDR_cis.R $OUT.txt.gz 0.05 $OUT
done
