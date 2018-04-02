#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools
QTLtools=/home/vitor/qtltools/bin/QTLtools

sample=$1

indexDIR=/home/vitor/hlaexpression/index_genome/star/index
fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN 8 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outSAMunmapped Within KeepPairs\
  --outSAMtype BAM SortedByCoordinate\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.sortedByCoord.out.bam
bai=$bam.bai
gencode=/home/vitor/gencode_data/gencode.v25.primary_assembly.annotation.gtf.gz
outQuant=./quantifications/$sample

$samtools index $bam $bai

$QTLtools quan --bam $bam --gtf $gencode --sample $sample \
    --out-prefix $outQuant \
    --filter-mapping-quality 255 --filter-mismatch 6 --filter-mismatch-total 6 \
    --tpm --check-proper-pairing --check-consistency --no-merge \
    --filter-remove-duplicates

rm ${outPrefix}*
