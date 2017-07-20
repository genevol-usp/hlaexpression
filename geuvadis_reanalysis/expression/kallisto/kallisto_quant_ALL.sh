#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq
fastqR1=$fastq_dir/$sample\_1.fastq.gz 
fastqR2=$fastq_dir/$sample\_2.fastq.gz 
index=./index/gencode.v26.ALL.transcripts.idx
outdir=./quantifications_ALL

$kallisto quant -i $index -t 1 -o $outdir/$sample --bias $fastqR1 $fastqR2\
  &> $log