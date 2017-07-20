#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastqR1=../../data/fastq/$sample\_1.fq.gz
fastqR2=../../data/fastq/$sample\_2.fq.gz

index=../../../geuvadis_reanalysis/expression/kallisto/index/gencode.v26.CHR.transcripts.idx

outdir=./quantifications_CHR
log=$outdir/log/$sample.quant.log

$kallisto quant -i $index -t 1 -o $outdir/$sample --bias $fastqR1 $fastqR2\
  &> $log 
