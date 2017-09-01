#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=./index
fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outMap=./mappings_1
outQuant=./quantifications_1
outPrefix=$outMap/$sample\_

mkdir -p $outMap
mkdir -p $outQuant

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 1\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 3000\
  --winAnchorMultimapNmax 6000\
  --alignIntronMax 1\
  --alignEndsType EndToEnd\
  --outSAMunmapped None\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=../../../imgt_index/gencode.v25.PRI.IMGT.transcripts.fa
out=$outQuant/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

rm $outPrefix*
