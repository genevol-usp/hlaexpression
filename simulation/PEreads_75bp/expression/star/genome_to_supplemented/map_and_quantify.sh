#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools
seqtk=/home/vitor/seqtk/seqtk
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/index_genome/star/index
fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_1
outPrefix=${outMap}/${sample}_

#$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
#  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
#  --outSAMunmapped Within KeepPairs\
#  --outSAMtype BAM SortedByCoordinate\
#  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.sortedByCoord.out.bam
hlabam=${outPrefix}hla.bam
unmapbam=${outPrefix}unmapped.bam
finalbam=${outPrefix}hla_and_unmapped.bam

#$samtools view -f 0x2 -F 0x100 -b -L hla.bed $bam > $hlabam
#$samtools view -F 0x2 $bam > $unmapbam
#$samtools merge $finalbam $hlabam $unmapbam

readid=${outPrefix}readids.txt
fq12=${outPrefix}1.fq
fq22=${outPrefix}2.fq

#$samtools view $finalbam | cut -f1 | sort | uniq > $readid
#
#$seqtk subseq $fq1 $readid > $fq12 
#$seqtk subseq $fq2 $readid > $fq22 

indexDIR2=/home/vitor/hlaexpression/index_imgtonly/index

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR2\
    --readFilesIn $fq12 $fq22\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 10000\
    --winAnchorMultimapNmax 20000\
    --alignIntronMax 0\
    --alignEndsType EndToEnd\
    --outSAMunmapped None\
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

#rm ${outPrefix}*
