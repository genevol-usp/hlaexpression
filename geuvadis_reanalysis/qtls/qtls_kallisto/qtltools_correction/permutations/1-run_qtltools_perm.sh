#!/bin/bash

QTLtools_dir=/home/vitor/QTLtools
QTLtools=$QTLtools_dir/QTLtools_1.1_Ubuntu16.04_x86_64
parallel=/home/vitor/parallel

CHUNKS=64
SAMPLES=../../../genotypes/samples.eur
VCF=../../../genotypes/eur_maf05.vcf.gz
COV=../../../pca_genotypes/covariates_genos.txt
CMD_FILE=./cmd.txt
  
for pc in $(seq 0 5 100)
do  
  BED=../phenotypes/phenotypes_eur_$pc.bed.gz
  OUT=./results/permutations_$pc
  LOG=./log/pc$pc.log

  for j in $(seq 1 $CHUNKS)
  do
    echo $QTLtools cis --vcf $VCF --bed $BED --cov $COV\
      --include-samples $SAMPLES --normal --chunk $j $CHUNKS --permute 1000\
      --out $OUT\_$j.txt --log $LOG
  done
done>$CMD_FILE
  
$parallel --gnu -j 64 :::: $CMD_FILE

rm $CMD_FILE

for pc in $(seq 0 5 100)
do
  OUT=./results/permutations_$pc
  cat $OUT\_*.txt | gzip -c > $OUT.txt.gz
  rm $OUT\_*.txt

  Rscript $QTLtools_dir/script/runFDR_cis.R $OUT.txt.gz 0.05 $OUT
done