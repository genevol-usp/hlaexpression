#!/bin/bash

QTLtools_dir=/home/vitor/QTLtools
QTLtools=$QTLtools_dir/QTLtools
parallel=/home/vitor/parallel

JOBS=60
SAMPLES=../genotypes/samples.eur
VCF=../genotypes/eur_maf05.vcf.gz
COV=../covariates_genos.txt
CMD_FILE=./cmd.txt

for k in $(seq 0 5 10; seq 20 5 30; echo 40; seq 55 5 100)
do
  BED=./phenotypes/phenotypes_eur_$k.bed.gz
  OUT=./permutations/permutations_$k
  LOG=./log/k$k.log

  for j in $(seq 1 $JOBS)
  do
    echo $QTLtools cis --vcf $VCF --bed $BED --cov $COV \
      --include-samples $SAMPLES --normal --chunk $j $JOBS --permute 1000 \
      --out $OUT\_$j.txt --log $LOG
  done
done>$CMD_FILE

$parallel --gnu -j $JOBS :::: $CMD_FILE

rm $CMD_FILE

for k in $(seq 0 5 10; seq 20 5 30; echo 40; seq 55 5 100)
do
  OUT=./permutations/permutations_$k
  cat $OUT\_*.txt | gzip -c > $OUT.txt.gz
  rm $OUT\_*.txt

  Rscript $QTLtools_dir/script/runFDR_cis.R $OUT.txt.gz 0.05 $OUT
done
