#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1

indexDIR=/home/vitor/hlaexpression/imgt_index/star/indexPRI
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=$outMap/${sample}_
outQuant=./quantifications

$STAR --runMode alignReads --runThreadN 4 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 150\
  --winAnchorMultimapNmax 300\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.PRI.transcripts.fa
out=$outQuant/$sample

if [ -d "$out" ]; then
    rmdir $out
fi

$salmon quant -t $fasta -l IU -a $bam -o $out -p 4 --seqBias --gcBias

header=${outPrefix}header.sam
sampledir=$outMap/$sample
imgtbam=$sampledir/imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 -F 0x100 $bam |\
    grep -F -f ../../data/ids_to_filter.txt - |\
    cat $header - |\
    $samtools view -Sb - > $imgtbam

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ* ${outPrefix}header.sam
