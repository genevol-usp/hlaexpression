#!/bin/bash

write_custom_index=/home/vitor/hlaexpression/imgt_index/write_genotyped_alleles.R

outdir=./sample_indices

if [ -d "$outdir" ]; then
    rm -r $outdir
fi
    
mkdir -p $outdir

Rscript $write_custom_index ./genotype_calls.tsv $outdir
