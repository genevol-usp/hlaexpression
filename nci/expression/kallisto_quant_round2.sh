#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
fastq_dir=/home/vitor/hlaexpression/nci/data/fastq
fastqR1=(`ls -v $fastq_dir/$sample*R1_001.fastq.gz`)
fastqR2=(`ls -v $fastq_dir/$sample*R2_001.fastq.gz`) 

gencode=../index/gencode.v26.CHR.transcripts.noIMGT.fa 
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa
sample_idx=./sample_indices/index_$sample.idx

outdir=./quantifications_2
log=$outdir/log/$sample.log

cat $gencode $sample_hla > $sample_fa 

$kallisto index -i $sample_idx $sample_fa &> $outdir/log/$sample.index.log

if [ "${#fastqR1[@]}" == 1 ]; then 
 $kallisto quant -i $sample_idx -t 8 -o $outdir/$sample --bias ${fastqR1[0]} ${fastqR2[0]} &>> $log
elif [ "${#fastqR1[@]}" == 2 ]; then 
 $kallisto quant -i $sample_idx -t 8 -o $outdir/$sample --bias ${fastqR1[0]} ${fastqR2[0]} ${fastqR1[1]} ${fastqR2[1]} &>> $log
else
  echo "wrong number of fastq files"
fi

rm $sample_idx $sample_fa
