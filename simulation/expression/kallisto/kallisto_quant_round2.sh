#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
fastqR1=../../data/fastq/$sample\_1.fq.gz
fastqR2=../../data/fastq/$sample\_2.fq.gz

gencode=../../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.transcripts.noIMGT.fa 
sample_hla=./index/sample_indices/hla_$sample.fa
sample_fa=./index/sample_indices/index_$sample.fa
sample_idx=./index/sample_indices/index_$sample.idx

outdir=./quantifications_2
sampledir=$outdir/$sample
bam=$sampledir/alignments.bam
log=$outdir/log/$sample.quant.log

mkdir -p $outdir/log
mkdir -p $sampledir

cat $gencode $sample_hla > $sample_fa

$kallisto index -i $sample_idx $sample_fa &> $outdir/log/$sample.index.log

$kallisto quant -i $sample_idx -t 1 -o $sampledir --bias --pseudobam\
  $fastqR1 $fastqR2 2> $log | $samtools view -Sb - > $bam

rm $sample_fa $sample_idx
