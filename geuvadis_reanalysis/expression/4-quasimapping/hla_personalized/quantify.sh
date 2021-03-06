#!/bin/bash

salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1
CPUS=6
indexDIR=./sample_indices/$sample
fasta=/home/vitor/hlaexpression/imgt_index_v2/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=../../2-hla_typing/sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$salmon index -t $sample_fa -i $indexDIR --type quasi -k 31

fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outQuant=./quantifications
out=$outQuant/$sample

mkdir -p $outQuant

$salmon quant -i $indexDIR -l IU -1 $fq1 -2 $fq2 -o $out -p $CPUS --seqBias --gcBias

rm -r $sample_fa $indexDIR
