#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastqR1=./fastq/$sample\_1.fq.gz
fastqR2=./fastq/$sample\_2.fq.gz

index=../index/gencode.v26.CHR.IMGT.transcripts.idx

outdir=./quantifications_1
log=$outdir/log/$sample.log

$kallisto quant -i $index -t 1 -o $outdir/$sample $fastqR1 $fastqR2 &> $log 
