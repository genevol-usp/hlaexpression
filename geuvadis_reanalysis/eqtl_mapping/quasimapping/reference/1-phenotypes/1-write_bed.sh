#!/bin/bash

BED=../../../../expression/4-quasimapping/reference/quantifications_expressed50%.bed
OUT=./phenotypes_eur.bed

cp $BED $OUT
bgzip $OUT && tabix -p bed $OUT.gz
