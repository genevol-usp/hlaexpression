#!/bin/bash

write_custom_index=/home/vitor/hlaexpression/imgt_index_utr/write_genotyped_alleles.R

outdir=./sample_indices

mkdir -p $outdir

Rscript $write_custom_index ./quantifications_1/processed_imgt_quants.tsv $outdir