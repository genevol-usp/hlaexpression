#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

$QTLtools pca --bed ./phenotypes_eur.bed.gz --out ./phenotypes_eur_pcs --center --scale
