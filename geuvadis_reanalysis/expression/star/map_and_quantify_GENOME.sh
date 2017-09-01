#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=../../../imgt_index/star/indexGENOME
fq1=../../data/fastq/$sample\_1.fastq.gz
fq2=../../data/fastq/$sample\_2.fastq.gz
outMap=./mappings_GENOME
outQuant=./quantifications_GENOME

mkdir -p $outMap
mkdir -p $outQuant

outPrefix=$outMap/$sample\_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 10\
  --winAnchorMultimapNmax 50\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM SortedByCoordinate\
  --outFileNamePrefix $outPrefix
