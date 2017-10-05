#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastqR1=../../data/fastq/${sample}_1.fastq.gz
fastqR2=../../data/fastq/${sample}_2.fastq.gz
index=../../../../imgt_index/kallisto/gencode.v25.PRI.IMGT.transcripts.idx
outdir=./quantifications_1
sampledir=$outdir/$sample

$kallisto quant -i $index -t 1 -o $sampledir $fastqR1 $fastqR2 

awk 'NR==1 || $1 ~ /IMGT/ {print $1"\t"$4"\t"$5}' $sampledir/abundance.tsv >\
    $sampledir/abundance_imgt.tsv