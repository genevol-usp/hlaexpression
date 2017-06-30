#!/bin/bash

#Rscript make_imgt_index.R
#Rscript make_gencode_fastas.R 

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto 
fasta=gencode.v26.CHR.IMGT.transcripts.fa
out=gencode.v26.CHR.IMGT.transcripts.idx

$kallisto index -i $out $fasta
