#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastqR1=../../data/fastq/${sample}_1.fastq.gz
fastqR2=../../data/fastq/${sample}_2.fastq.gz

gencode=../../../imgt_index/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa
sample_idx=./sample_indices/index_$sample.idx

outdir=quantifications_2
sampledir=$outdir/$sample

cat $gencode $sample_hla > $sample_fa

$kallisto index -i $sample_idx $sample_fa

$kallisto quant -i $sample_idx -t 1 -o $sampledir --bias $fastqR1 $fastqR2

rm $sample_fa $sample_idx
