#!/bin/bash

input_fasta=$1
output_fastq=$(echo ${input_fasta/fasta/fq})

cat $input_fasta | \
  awk -v mate=$2 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {
		    print \
		    "@"gensub(/;/, " ", 1, gensub(/\//, "_", 1, $1))"/"mate"\n"\
		    $2"\n"\
		    "+\n"\
		    gensub(/./, ">", "g", $2)\
		  }' | gzip -c > $output_fastq.gz
		
rm $input_fasta
