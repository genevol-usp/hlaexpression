#!/bin/bash

STAR=/home/vitor/STAR

indexDIR=./indexCHR
fasta=/home/vitor/gencode_data/gencode.v25.transcripts.fa

mkdir -p $indexDIR

$STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $fasta --genomeChrBinNbits 10 --genomeSAindexNbases 13
