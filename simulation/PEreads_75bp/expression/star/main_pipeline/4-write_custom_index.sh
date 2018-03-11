#!/bin/bash

write_custom_index=/home/vitor/hlaexpression/imgt_index_v2/write_genotyped_alleles.R

outdir=./sample_indices

mkdir -p $outdir

Rscript $write_custom_index ./quantifications_MHC/processed_imgt_quants.tsv $outdir
