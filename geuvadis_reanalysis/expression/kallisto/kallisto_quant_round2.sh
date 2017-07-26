#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq
fastqR1=$fastq_dir/$sample\_1.fastq.gz
fastqR2=$fastq_dir/$sample\_2.fastq.gz

gencode=./index/gencode.v25.CHR.transcripts.noIMGT.fa  
sample_hla=./index/sample_indices/hla_$sample.fa
sample_fa=./index/sample_indices/index_$sample.fa
sample_idx=./index/sample_indices/index_$sample.idx

outdir=quantifications_2
sampledir=$outdir/$sample

cat $gencode $sample_hla > $sample_fa

$kallisto index -i $sample_idx $sample_fa &> $outdir/log/$sample.index.log

$kallisto quant -i $sample_idx -t 1 -o $sampledir --bias $fastqR1 $fastqR2\
  2> $outdir/log/$sample.quant.log

rm $sample_fa $sample_idx
