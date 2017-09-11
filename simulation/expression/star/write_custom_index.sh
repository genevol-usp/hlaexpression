#!/bin/bash

write_custom_index=../../../imgt_index/write_genotyped_alleles.R

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./sample_indices
