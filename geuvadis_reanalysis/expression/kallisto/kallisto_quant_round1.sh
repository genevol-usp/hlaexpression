#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq
fastqR1=$fastq_dir/$sample\_1.fastq.gz
fastqR2=$fastq_dir/$sample\_2.fastq.gz

index=./index/gencode.v25.CHR.IMGT.transcripts.idx

outdir=./quantifications_1
log=$outdir/log/$sample.quant.log

$kallisto quant -i $index -t 1 -o $outdir/$sample $fastqR1 $fastqR2 &> $log
