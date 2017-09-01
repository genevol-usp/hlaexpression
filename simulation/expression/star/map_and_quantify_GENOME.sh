#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=../../../imgt_index/star/indexGENOME
fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outMap=./mappings_GENOME
outQuant=./quantifications_GENOME
outPrefix=$outMap/$sample\_

mkdir -p $outMap
mkdir -p $outQuant

$STAR --runMode alignReads --runThreadN 16 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.06\
  --outFilterMultimapScoreRange 0\
  --outFilterType BySJout\
  --outFilterMultimapNmax 20\
  --winAnchorMultimapNmax 50\
  --alignIntronMax 0\
  --alignEndsType Local\
  --alignSJoverhangMin 8\
  --alignSJDBoverhangMin 3\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM SortedByCoordinate\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix
