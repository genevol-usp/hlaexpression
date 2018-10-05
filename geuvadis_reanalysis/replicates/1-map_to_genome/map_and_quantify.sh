#!/bin/bash

STAR=/home/vitor/Libraries/STAR
samtools=/home/vitor/Libraries/samtools-1.3.1/samtools
seqtk=/home/vitor/Libraries/seqtk/seqtk
salmon=/home/vitor/Libraries/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=8
indexDIR=/home/vitor/hlaexpression/index_genome/star/index
outMap=./mappings
runid=$(pwgen 5 1)
outPrefix=/scratch/genevol/users/vitor/${runid}_${sample}_

fq1=/scratch/genevol/users/vitor/geuvadis_replicates/${sample}_1.fastq.gz
fq2=/scratch/genevol/users/vitor/geuvadis_replicates/${sample}_2.fastq.gz

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 20\
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within KeepPairs\
  --outFileNamePrefix $outPrefix

bamGenome=${outPrefix}Aligned.sortedByCoord.out.bam
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt

$samtools index $bamGenome $bamGenome.bai

$samtools view $bamGenome "chr6:29722000-33144000" |\
    cut -f1 |\
    sort |\
    uniq > $readsalign

$samtools view -F 0x2 $bamGenome |\
    cut -f1 |\
    sort |\
    uniq > $readsunmap

cat $readsalign $readsunmap |\
    sort |\
    uniq > $readids

fq12=./mappings/mhc_fqs/${sample}_1.fq
fq22=./mappings/mhc_fqs/${sample}_2.fq

$seqtk subseq $fq1 $readids > $fq12 
$seqtk subseq $fq2 $readids > $fq22 

rm ${outPrefix}*
