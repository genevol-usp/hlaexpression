#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools

$QTLtools pca --bed ./phenotypes_eur.bed.gz --out ./phenotypes_eur_pcs\
  --center --scale
