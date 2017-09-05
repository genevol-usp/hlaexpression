#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1

indexDIR=../../../imgt_index/star/indexGENOME
fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outMap=./mappings_GENOME
outQuant=./quantifications_GENOME
outPrefix=$outMap/$sample\_

mkdir -p $outMap
mkdir -p $outQuant

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterType BySJout\
  --outFilterMultimapNmax 20\
  --winAnchorMultimapNmax 50\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype None\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.toTranscriptome.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.PRI.transcripts.fa
out=$outQuant/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

header=${outPrefix}header.sam
imgtbam=${outPrefix}imgt.bam

$samtools view -H $bam > $header

$samtools view $bam |\
  LC_ALL=C grep -F -f ./ids_to_filter.txt |\
  cat $header - |\
  $samtools view -Sb - |\
  $samtools view -b -f 0x2 -F 0x100 - > $imgtbam

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ* ${outPrefix}header.sam
