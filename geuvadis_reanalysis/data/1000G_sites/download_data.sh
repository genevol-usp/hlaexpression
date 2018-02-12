#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools

#for i in {1..22}
#do
#	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${i}_GRCh38_sites.20170504.vcf.gz &
#done

$bcftools concat -o all_autosomes_sites.vcf.gz -O z $(ls -v ALL.chr*.vcf.gz) 
