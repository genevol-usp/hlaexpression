#!/bin/bash

PC=40
OUT=./conditional_$PC

cat ${OUT}_*.txt | gzip -c > ${OUT}_all.txt.gz
rm ${OUT}_*.txt 
