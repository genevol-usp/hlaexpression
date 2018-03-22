#!/bin/bash

caveman=/home/vitor/CaVEMaN
EQTL=./eqtl.list
BED=./phenotypes_with_eqtl.bed
VCF=/home/vitor/hlaexpression/geuvadis_reanalysis/qtls/genotypes/eur_maf05.vcf.gz
COV=/home/vitor/hlaexpression/geuvadis_reanalysis/qtls/pca_genotypes/covariates_genos.txt
OUT=./phenotypes_corrected.bed

$caveman --single-signal --eqtl $EQTL --bed $BED --vcf $VCF --cov $COV --out $OUT
