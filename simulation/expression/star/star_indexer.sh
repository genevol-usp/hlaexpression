#!/bin/bash

STAR=/home/vitor/STAR

indexDIR=./index
fasta=../../../geuvadis_reanalysis/expression/kallisto/index/gencode.v25.CHR.IMGT.transcripts.fa

$STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $fasta --genomeChrBinNbits 10 --genomeSAindexNbases 13
