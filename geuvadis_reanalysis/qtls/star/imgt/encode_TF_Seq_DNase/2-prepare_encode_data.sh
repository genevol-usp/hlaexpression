#!/bin/bash

liftDIR=/home/vitor/liftOver
liftover=$liftDIR/liftOver
chain=$liftDIR/hg19ToHg38.over.chain

# Files obtained from O. Delaneau:
TF=/home/vitor/ENCODE/TranscriptionFactors/step5_merged/ALL.bed.gz
SEG=/home/vitor/ENCODE/Segmentation/step0_raw/wgEncodeAwgSegmentationCombinedGm12878.bed
DNase=/home/vitor/ENCODE/DNAse/step1_bed/wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.bed.gz
HM=/home/vitor/ENCODE/HistoneModifications/step2_merge/HM.ALL.bed.gz

tmp1=TF.ENCODE.chr6.bed
tmp2=Seg.ENCODE.chr6.bed 
tmp3=DNase.ENCODE.chr6.bed
tmp4=HM.ENCODE.chr6.bed

out1=$( echo $tmp1 | sed 's/bed/hg38.bed/' )
out2=$( echo $tmp2 | sed 's/bed/hg38.bed/' )
out3=$( echo $tmp3 | sed 's/bed/hg38.bed/' )
out4=$( echo $tmp4 | sed 's/bed/hg38.bed/' )

zcat $TF | awk '$1 == "chr6"' > $tmp1
cat $SEG | awk '$1 == "chr6"' > $tmp2 
zcat $DNase | awk '$1 == "chr6"' > $tmp3
zcat $HM | awk '$1 == "chr6"' > $tmp4

$liftover $tmp1 $chain $out1 ./unmapped_TF
$liftover $tmp2 $chain $out2 ./unmapped_Seg
$liftover $tmp3 $chain $out3 ./unmapped_DNase
$liftover $tmp4 $chain $out4 ./unmapped_DNase

rm $tmp1 $tmp2 $tmp3 $tmp4
