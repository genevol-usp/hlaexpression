#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools
seqtk=/home/vitor/seqtk/seqtk
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/index_genome/star/index
fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outSAMtype None\
  --outReadsUnmapped Fastx\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.toTranscriptome.out.bam
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt

$samtools view -L imgt.bed $bam |\
    cut -f1 |\
    sort |\
    uniq > $readsalign

sed -n '1~4p' ${outPrefix}Unmapped.out.mate1 |\
    cut -f1 |\
    sed 's/^@//g' |\
    sort |\
    uniq > $readsunmap

cat $readsalign $readsunmap |\
    sort |\
    uniq > $readids

fq12=${outPrefix}1.fq
fq22=${outPrefix}2.fq

$seqtk subseq $fq1 $readids > $fq12 
$seqtk subseq $fq2 $readids > $fq22 

indexDIR2=/home/vitor/hlaexpression/index_imgtonly/index

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR2\
    --readFilesIn $fq12 $fq22\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 5000\
    --winAnchorMultimapNmax 10000\
    --alignIntronMax 0\
    --alignEndsType EndToEnd\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}step2_

fasta=/home/vitor/hlaexpression/imgt_index_v2/imgt_index.fa
bam2=${outPrefix}step2_Aligned.out.bam
out=./quantifications_1/$sample

if [ -d "$out" ]; then
    rm -r $out
fi

$salmon quant -t $fasta -l IU -a $bam2 -o $out -p 6 --seqBias --gcBias

rm ${outPrefix}*
