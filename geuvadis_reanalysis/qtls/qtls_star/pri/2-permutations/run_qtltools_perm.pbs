#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -l walltime=24:00:00 
#PBS -q short 
#PBS -j oe 
#PBS -o /home/vitor/hlaexpression/geuvadis_reanalysis/qtls/qtls_star/imgt/2-permutations/log/$PBS_JOBID.log

cd $PBS_O_WORKDIR

QTLtools_dir=/home/vitor/QTLtools
QTLtools=$QTLtools_dir/QTLtools_1.1_Ubuntu16.04_x86_64
parallel=/home/vitor/parallel

SAMPLES=../../genotypes/samples.eur
VCF=../../genotypes/eur_maf05.vcf.gz
COV=../../pca_genotypes/covariates_genos.txt
CMD_FILE=./cmd.txt
  
BED=../phenotypes/phenotypes_eur_$PC.bed.gz
OUT=./results/permutations_$PC
LOG=./log/pc$PC.log

$QTLtools cis --vcf $VCF --bed $BED --cov $COV\
  --include-samples $SAMPLES --normal --chunk $CHUCK 60 --permute 1000\
  --out ${OUT}_$CHUNK.txt
  
