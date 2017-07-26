#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=../data/fastq
fastqR1=(`ls -v $fastq_dir/$sample*R1_001.fastq.gz`)
fastqR2=(`ls -v $fastq_dir/$sample*R2_001.fastq.gz`)

index=../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.IMGT.transcripts.idx

outdir=./quantifications_1
log=$outdir/log/$sample.log

mkdir -p $outdir/log

if [ "${#fastqR1[@]}" == 1 ]; then
  $kallisto quant -i $index -t 1 -o $outdir/$sample ${fastqR1[0]} ${fastqR2[0]} &> $log
elif [ "${#fastqR1[@]}" == 2 ]; then
  $kallisto quant -i $index -t 1 -o $outdir/$sample ${fastqR1[0]} ${fastqR2[0]} ${fastqR1[1]} ${fastqR2[1]} &> $log
else
  echo "wrong number of fastq files"
fi
