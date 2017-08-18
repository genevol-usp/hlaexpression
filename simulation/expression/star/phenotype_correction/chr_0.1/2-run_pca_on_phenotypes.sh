#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

$QTLtools pca --bed ./phenotypes.bed.gz --out ./phenotypes_pcs --center --scale
