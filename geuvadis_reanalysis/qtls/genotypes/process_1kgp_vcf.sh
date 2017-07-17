#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf_dir=/home/vitor/hlaexpression/geuvadis_reanalysis/data/1000G_genotypes
samples=./samples.all 
out=./genotypes_maf05_either_eur_yri.vcf.gz

for i in $(ls $vcf_dir/ALL.chr*gz) 
do 
  $bcftools view --samples-file $samples\
    -o $(echo $i | sed "s|$vcf_dir|.|" | sed 's|.vcf.gz|.GeuvSamples.vcf.gz|')\
    -O z $i & 
done
wait

for i in $(ls ALL.chr*GeuvSamples.vcf.gz)
do 
  $bcftools view --genotype ^miss -m2 -M2\
    -o $(echo $i | sed 's|.GeuvSamples.vcf.gz|.GeuvSamples_nomiss_bi.vcf.gz|')\
    -O z $i &
done
wait

$bcftools concat -o concat.vcf.gz -O z\
  $(ls -v ALL.chr*GeuvSamples_nomiss_bi.vcf.gz)

zcat concat.vcf.gz | grep -v "^#" | awk '{print $3}' | uniq -c |\
  awk '$1>=2{print $2}' > duplicate_ids.txt

$bcftools view --exclude '%ID=@duplicate_ids.txt'\
  -o concat_nodups.vcf.gz -O z concat.vcf.gz

$bcftools view --samples-file samples.eur --min-af 0.05:minor\
  -o eur_maf05.vcf.gz -O z concat_nodups.vcf.gz

$bcftools view --samples-file samples.yri --min-af 0.05:minor\
  -o yri_maf05.vcf.gz -O z concat_nodups.vcf.gz

zcat eur_maf05.vcf.gz yri_maf05.vcf.gz | grep -v "^#" | cut -f3 | sort |\
  uniq > final_variants.txt

$bcftools view --include '%ID=@final_variants.txt'\
  -o $out -O z concat_nodups.vcf.gz

rm *GeuvSamples_nomiss_bi.vcf.gz *GeuvSamples.vcf.gz concat.vcf.gz\
  concat_nodups.vcf.gz

tabix -p vcf $out
tabix -p vcf eur_maf05.vcf.gz
