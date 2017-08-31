#!/bin/bash

salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1
index=./index
fastq1=../../data/fastq/$sample\_1.fq.gz
fastq2=../../data/fastq/$sample\_2.fq.gz
outQuant=./quantifications_1
out=$outQuant/$sample

mkdir -p $outQuant

$salmon quant -i $index -l IU -1 $fastq1 -2 $fastq2 -o $out -p 6 \
  --maxReadOcc 4000 --useVBOpt
  #--seqBias --gcBias --posBias --maxReadOcc 4000
