#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

transcripts_chr=gencode.v26.CHR.transcripts
transcripts_all=gencode.v26.ALL.transcripts

kallisto index -i $transcripts_chr.idx $transcripts_chr.fa
kallisto index -i $transcripts_all.idx $transcripts_all.fa
