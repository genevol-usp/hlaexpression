#!/bin/bash

write_custom_index=/home/vitor/hlaexpression/imgt_index_v2/write_genotyped_alleles.R

outdir=./sample_indices

if [ -d "$outdir" ]; then
    rm -r $outdir
    mkdir -p $outdir
fi

Rscript $write_custom_index ./quantifications_genos/processed_imgt_quants.tsv $outdir
