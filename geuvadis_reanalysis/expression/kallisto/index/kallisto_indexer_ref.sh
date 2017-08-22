#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

transcripts_chr=gencode.v25.CHR.transcripts
transcripts_all=gencode.v25.ALL.transcripts

$kallisto index -i $transcripts_chr.idx ../../../../imgt_index/$transcripts_chr.fa&
$kallisto index -i $transcripts_all.idx ../../../../imgt_index/$transcripts_all.fa
