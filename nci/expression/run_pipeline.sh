#!/bin/bash

parallel=/home/vitor/parallel

write_custom_index=../../geuvadis_reanalysis/expression/kallisto/index/write_genotyped_alleles.R

$parallel --gnu -j 16 ./kallisto_quant_round1.sh {} :::: nci_samples.txt

Rscript process_kallisto.R 1

Rscript $write_custom_index ./quantifications_1/processed_quant.tsv ./index/sample_indices

$parallel --gnu -j 64 ./kallisto_quant_round2_bam.sh {} :::: ./nci_samples.txt

Rscript process_kallisto.R 2
