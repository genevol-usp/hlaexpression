#!/bin/bash

write_custom_index=../../../imgt_index/write_genotyped_alleles.R

out=./sample_indices
mkdir -p $out

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv $out
