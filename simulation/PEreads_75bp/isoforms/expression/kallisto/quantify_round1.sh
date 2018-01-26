#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fq1=../../tmp_data/fastq/${sample}_1.fastq.gz
fq2=../../tmp_data/fastq/${sample}_2.fastq.gz
index=/home/vitor/hlaexpression/imgt_index/kallisto/gencode.v25.PRI.IMGT.transcripts.idx
outdir=./quantifications_1
sampledir=$outdir/$sample

$kallisto quant -i $index -t 1 -o $sampledir $fq1 $fq2
