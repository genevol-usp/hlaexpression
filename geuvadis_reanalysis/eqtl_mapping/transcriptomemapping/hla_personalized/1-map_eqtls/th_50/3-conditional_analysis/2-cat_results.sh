#!/bin/bash

PC=60
OUT=./conditional_$PC

cat ${OUT}_*.txt | gzip -c > ${OUT}_all.txt.gz
rm ${OUT}_*.txt 
