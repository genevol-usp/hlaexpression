#!/bin/bash

parallel=/home/vitor/parallel

#mkdir -p ./fastq
#mkdir -p ./fasta_{01..50}

#$parallel --gnu -j 50 Rscript 2-polyester.R {} ::: {1..50} 

#for i in {1..2}
#do
#  for j in {01..50}
#  do
#    cat ./fasta_$j/sample_01_$i.fasta |\
#    awk -v mate=$i 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {
#      print\
#      "@"gensub(/;/, " ", 1, gensub(/\//, "_", 1, $1))"/"mate"\n"\
#      $2"\n"\
#      "+\n"\
#      gensub(/./, ">", "g", $2)\
#    }' | gzip -c > ./fastq/sample_$j\_$i.fq.gz &
#  done
#  wait
#done
#wait
#
rm -r ./fasta_{01..50}
