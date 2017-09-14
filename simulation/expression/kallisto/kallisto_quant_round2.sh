#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
fastqR1=../../data/fastq/$sample\_1.fastq.gz
fastqR2=../../data/fastq/$sample\_2.fastq.gz

gencode=../../../imgt_index/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa
sample_idx=./sample_indices/index_$sample.idx

outdir=./quantifications_2
sampledir=$outdir/$sample
bam=$sampledir/alignments.bam

mkdir -p $sampledir

cat $gencode $sample_hla > $sample_fa

$kallisto index -i $sample_idx $sample_fa

$kallisto quant -i $sample_idx -t 1 -o $sampledir --bias --pseudobam\
  $fastqR1 $fastqR2 | $samtools view -Sb - > $bam

header=$sampledir/header.sam
imgtbam=$sampledir/imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 $bam |\
  awk -F $'\t' '$1 ~ /IMGT/ || $3 ~ /IMGT/' |\
  cat $header - |\
  $samtools view -Sb - > $imgtbam

rm $sample_fa $sample_idx $bam
