#!/bin/bash

write_custom_index=/home/vitor/hlaexpression/imgt_index_v2/write_genotyped_alleles.R

outdir=./sample_indices

if [ -d "$outdir" ]; then
    rm -r $outdir
    mkdir -p $outdir
fi

mkdir -p $outdir

Rscript $write_custom_index ./quantifications_MHC/imgt_quants_top5.tsv $outdir
