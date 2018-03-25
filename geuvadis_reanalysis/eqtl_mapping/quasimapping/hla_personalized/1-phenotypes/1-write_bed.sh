#!/bin/bash

BED=/home/vitor/hlaexpression/geuvadis_reanalysis/expression/star/main_pipeline/quantifications_expressed50%.bed
OUT=./phenotypes_eur.bed

cp $BED $OUT
bgzip $OUT && tabix -p bed $OUT.gz
