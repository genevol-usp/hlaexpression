#!/bin/bash

salmon=/home/vitor/Libraries/Salmon-latest_linux_x86_64/bin/salmon

sample=$1
CPUS=8
indexDIR=/home/vitor/hlaexpression/index_transcriptome/salmon/index

fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outQuant=./quantifications
out=$outQuant/$sample

mkdir -p $outQuant

$salmon quant -i $indexDIR -l IU -1 $fq1 -2 $fq2 -o $out -p $CPUS --seqBias --gcBias
