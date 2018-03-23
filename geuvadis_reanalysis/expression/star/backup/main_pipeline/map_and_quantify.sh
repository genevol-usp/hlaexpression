#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools
seqtk=/home/vitor/seqtk/seqtk
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=8
indexDIR=/home/vitor/hlaexpression/index_genome/star/index
fq1=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_1.fastq.gz
fq2=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 20\
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within KeepPairs\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bamTransc=${outPrefix}Aligned.toTranscriptome.out.bam
fastaTransc=/home/vitor/hlaexpression/index_transcriptome/gencode.v25.PRI.transcripts.fa
outTransc=./quantifications_transcriptome/$sample

if [ -d "$outTransc" ]; then
    rm -r $outTransc
fi

$salmon quant -t $fastaTransc -l IU -a $bamTransc -o $outTransc -p $CPUS --seqBias --gcBias

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

fq12=./mappings/fqs/${sample}_1.fq
fq22=./mappings/fqs/${sample}_2.fq

$seqtk subseq $fq1 $readids > $fq12 
$seqtk subseq $fq2 $readids > $fq22 

indexDIR2=/home/vitor/hlaexpression/index_mhc/star/index

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR2\
    --readFilesIn $fq12 $fq22\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}MHC_

fastaMHC=/home/vitor/hlaexpression/index_mhc/gencode.v25.MHC.IMGT.transcripts.fa
bamMHC=${outPrefix}MHC_Aligned.out.bam
outMHC=./quantifications_MHC/$sample

if [ -d "$outMHC" ]; then
    rm -r $outMHC
fi

$salmon quant -t $fastaMHC -l IU -a $bamMHC -o $outMHC -p $CPUS

rm ${outPrefix}*
