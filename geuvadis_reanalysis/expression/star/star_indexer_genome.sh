#!/bin/bash

STAR=/home/vitor/STAR

indexDIR=./indexGENOME
genome=/home/vitor/gencode_data/GRCh38.p7.genome.fa
gtf=/home/vitor/gencode_data/gencode.v25.annotation.gtf
overhang=75

mkdir -p $indexDIR

zcat $genome.gz > $genome
zcat $gtf.gz > $gtf

$STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbOverhang $overhang 

rm $genome $gtf
