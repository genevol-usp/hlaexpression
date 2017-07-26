#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto 
fasta=gencode.v25.CHR.IMGT.transcripts.fa
out=gencode.v25.CHR.IMGT.transcripts.idx

$kallisto index -i $out $fasta
