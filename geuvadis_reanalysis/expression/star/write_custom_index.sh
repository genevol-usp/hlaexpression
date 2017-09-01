#!/bin/bash

write_custom_index=../../../geuvadis_reanalysis/expression/kallisto/index/write_genotyped_alleles.R

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./sample_indices
