#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=/home/vitor/hlaexpression/nci/data/fastq
fastqR1=(`ls -v $fastq_dir/$sample*R1_001.fastq.gz`)
fastqR2=(`ls -v $fastq_dir/$sample*R2_001.fastq.gz`)

index=../index/gencode.v26.CHR.IMGT.transcripts.idx

outdir=./quantifications_1
log=$outdir/log/$sample.log

if [ "${#fastqR1[@]}" == 1 ]; then
  $kallisto quant -i $index -t 8 -o $outdir/$sample ${fastqR1[0]} ${fastqR2[0]} &> $log
elif [ "${#fastqR1[@]}" == 2 ]; then
  $kallisto quant -i $index -t 8 -o $outdir/$sample ${fastqR1[0]} ${fastqR2[0]} ${fastqR1[1]} ${fastqR2[1]} &> $log
else
  echo "wrong number of fastq files"
fi