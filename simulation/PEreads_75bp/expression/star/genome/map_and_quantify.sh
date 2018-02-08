#!/bin/bash

STAR=/home/vitor/STAR

sample=$1

indexDIR=/home/vitor/hlaexpression/index_genome/index
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN 8 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outSAMtype BAM\
  --outFileNamePrefix $outPrefix
