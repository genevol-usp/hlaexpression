#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastqR1=../../data/fastq/${sample}_1.fastq.gz
fastqR2=../../data/fastq/${sample}_2.fastq.gz
index=../../../imgt_index/kallisto/gencode.v25.PRI.IMGT.transcripts.idx
outdir=./quantifications_1

$kallisto quant -i $index -t 1 -o $outdir/$sample $fastqR1 $fastqR2

