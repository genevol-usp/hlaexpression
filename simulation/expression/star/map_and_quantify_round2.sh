#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1
indexDIR=./sample_indices/$sample

mkdir -p $indexDIR

fasta=../../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.transcripts.noIMGT.fa
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $sample_fa --genomeChrBinNbits 11 --genomeSAindexNbases 13

fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outPrefix=./mappings_2/$sample\_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 50\
  --winAnchorMultimapNmax 100\
  --alignIntronMax 1\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

rm -r $indexDIR

bam=./mappings_2/$sample\_Aligned.out.bam
out=./quantifications_2/$sample

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p 6

rm $sample_fa $sample_hla
