#!/bin/bash

STAR=/home/vitor/STAR

indexDIR=./index
fasta=../kallisto/index/gencode.v26.CHR.IMGT.transcripts.fa

$STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $fasta --genomeChrBinNbits 10 --genomeSAindexNbases 13
