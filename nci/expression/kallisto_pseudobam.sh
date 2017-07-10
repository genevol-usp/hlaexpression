#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
fastqdir=/home/vitor/hlaexpression/nci/data/fastq
fastqR1=(`ls -v $fastqdir/$sample*R1_001.fastq.gz`)
fastqR2=(`ls -v $fastqdir/$sample*R2_001.fastq.gz`)

outdir=./alignments
sampledir=$outdir/$sample

mkdir -p ./alignments/log
mkdir -p $sampledir

gencode=../index/gencode.v26.CHR.transcripts.noIMGT.fa
hla=./sample_indices/hla_$sample.fa
index_fa=./sample_indices/index_$sample.fa
index_idx=./sample_indices/index_$sample.idx

log=$outdir/log/$sample.quant.log
outbam=$sampledir/alignments.bam

imgtbam=$sampledir/imgt.bam
imgtmapped=$sampledir/imgt_mapped.bam
imgtsorted=$sampledir/imgt_mapped_sorted.bam
bed=./imgt.bed
cov=./alignments/$sample/imgt.coverage

cat $gencode $hla > $index_fa

$kallisto index -i $index_idx $index_fa &> $outdir/log/$sample.index.log

if [ "${#fastqR1[@]}" == 1 ]; then
  $kallisto quant -i $index_idx -t 1 -o $sampledir --bias --pseudobam\
    ${fastqR1[0]} ${fastqR2[0]} 2> $log | $samtools view -Sb - > $outbam
elif [ "${#fastqR1[@]}" == 2 ]; then
  $kallisto quant -i $index_idx -t 1 -o $sampledir --bias --pseudobam\
    ${fastqR1[0]} ${fastqR2[0]} ${fastqR1[1]} ${fastqR2[1]} 2> $log |\
    $samtools view -Sb - > $outbam
else
  echo "wrong number of fastq files"
fi

$samtools view -H $outbam > $sampledir/header.sam

$samtools view $outbam | grep IMGT | cat $sampledir/header.sam - |\
  $samtools view -Sb - > $imgtbam

$samtools view -b -f 0x2 $imgtbam > $imgtmapped

$samtools sort $imgtmapped -o $imgtsorted

$samtools depth -a -m 100000 -b $bed $imgtsorted > $cov

rm $outbam $sampledir/header.sam $index_fa $index_idx $sampledir/abundance.* \
  $sampledir/run_info.json 
rm $imgtbam $imgtmapped
