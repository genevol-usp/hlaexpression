#!/bin/bash

#wget http://www.regulomedb.org/downloads/RegulomeDB.dbSNP141.txt.gz

zcat RegulomeDB.dbSNP141.txt.gz |\
    awk '$1 == "chr6" && $2 > 25000000 && X2 < 35000000' > RegulomeDB.dbSNP141.chr6.txt
