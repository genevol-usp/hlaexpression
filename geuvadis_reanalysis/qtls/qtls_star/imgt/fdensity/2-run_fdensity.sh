#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

qtl=./hla_qtls.bed
tfbs=../tfbs_encode/chr6_tfbs_hg38.bed
out=./density_tf_around_qtl.txt

$QTLtools fdensity --qtl $qtl --bed $tfbs --out $out
