#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
indexDIR=./sample_indices/$sample
fasta=../../../imgt_index/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $indexDIR\
    --genomeFastaFiles $sample_fa\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${indexDIR}_

fq1=../../data/fastq/$sample\_1.fastq.gz
fq2=../../data/fastq/$sample\_2.fastq.gz
outMap=./mappings_2
outQuant=./quantifications_2
outPrefix=$outMap/$sample\_

$STAR --runMode alignReads --runThreadN 16 --genomeDir $indexDIR\
    --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
    --outFilterMismatchNmax 999\
    --outFilterMismatchNoverReadLmax 0.04\
    --outFilterMultimapScoreRange 1\
    --outFilterMultimapNmax 100\
    --winAnchorMultimapNmax 200\
    --alignIntronMax 0\
    --alignEndsType Local\
    --outSAMunmapped Within KeepPairs\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix $outPrefix

rm -r $indexDIR ${indexDIR}_Log.out

bam=${outPrefix}Aligned.out.bam
out=$outQuant/$sample

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p 6

rm $sample_fa $sample_hla

header=${outPrefix}header.sam
imgtbam=${outPrefix}imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 -F 0x100 $bam |\
    awk -F $'\t' '$1 ~ /IMGT/ || $3 ~ /IMGT/' |\
    cat $header - |\
    $samtools view -Sb - > $imgtbam

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ* $header 
