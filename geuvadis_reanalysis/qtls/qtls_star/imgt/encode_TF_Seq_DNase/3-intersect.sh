#!/bin/bash

hla=./hla.qtls.bed
TF=./TF.ENCODE.chr6.hg38.bed
Seg=./Seg.ENCODE.chr6.hg38.bed
DNase=./DNase.ENCODE.chr6.hg38.bed

bedtools intersect -a $hla -b $TF > hla.TF.intersect.bed
bedtools intersect -b $hla -a $TF > TF.hla.intersect.bed

bedtools intersect -a $hla -b $Seg > hla.Seg.intersect.bed
bedtools intersect -b $hla -a $Seg > Seg.hla.intersect.bed

bedtools intersect -a $hla -b $DNase > hla.DNase.intersect.bed
bedtools intersect -b $hla -a $DNase > DNase.hla.intersect.bed

