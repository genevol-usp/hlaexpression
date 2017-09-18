#!/bin/bash

write_custom_index=../../../imgt_index/write_genotyped_alleles.R

mkdir -p ./sample_indices

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./sample_indices
