#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools

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

hlabam=${outPrefix}hla.bam
unmapbam=${outPrefix}unmapped.bam
finalbam=${outPrefix}hla_and_unmapped.bam

$samtools view -f 0x2 -F 0x100 -b -L hla.bed $bam > $hlabam
$samtools view -F 0x2 $bam > $unmapbam
$samtools merge $finalbam $hlabam $unmapbam


#bam=${outPrefix}Aligned.sortedByCoord.out.bam
#bai=$bam.bai
#gencode=/home/vitor/gencode_data/gencode.v25.primary_assembly.annotation.gtf.gz
#outQuant=./quantifications/$sample
#
#$samtools index $bam $bai
#
#$QTLtools quan --bam $bam --gtf $gencode --samples $sample \
#    --out-prefix $outQuant \
#    --filter-mapping-quality 255 --filter-mismatch 6 --filter-mismatch-total 6 \
#    --rpkm --check-proper-pairing --check-consistency --no-merge \
#    --filter-remove-duplicates
#
#header=${outPrefix}header.sam
#bamReadsGenerated=${outPrefix}ReadsGenerated.bam
#bamReadsMappedHLA=${outPrefix}ReadsMappedHLA.bam
#
#$samtools view -H $bam > $header
#
#$samtools view -f 0x2 -F 0x100 $bam |\
#    grep -F -f ../../../data/ids_to_filter.txt - |\
#    cat $header - |\
#    $samtools view -Sb - > $bamReadsGenerated
#
#$samtools view -b -L hla.bed $bam > $bamReadsMappedHLA
#
#rm $bam $bai $header ${outPrefix}Log.out ${outPrefix}Log.progress.out\
#    ${outPrefix}SJ*
