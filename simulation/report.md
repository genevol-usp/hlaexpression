Report
================

Notes:
======

Index nomenclature:

-   pri: Primary assembly of the reference genome, containing the reference chromosomes and scaffolds (no alternate haplotypes)

-   imgt: "pri" index supplemented with IMGT references

kallisto
========

Genotyping
----------

| locus    |  accuracy (%)|
|:---------|-------------:|
| HLA-A    |           100|
| HLA-B    |           100|
| HLA-C    |           100|
| HLA-DPB1 |            91|
| HLA-DQA1 |           100|
| HLA-DQB1 |           100|
| HLA-DRB1 |           100|

Expression
----------

<img src="./plots/kallisto_prop_mapped.png" width="1200" />

STAR + Salmon
=============

Genotyping
----------

| locus    |  accuracy (%)|
|:---------|-------------:|
| HLA-A    |           100|
| HLA-B    |           100|
| HLA-C    |           100|
| HLA-DPB1 |            92|
| HLA-DQA1 |           100|
| HLA-DQB1 |           100|
| HLA-DRB1 |           100|

Expression
----------

<img src="./plots/star_prop_mapped.png" width="1200" />

Quality assessment
==================

Percentage of simulated reads not aligned, aligned exclusively to a different reference, or aligned to original reference but as a multimap
-------------------------------------------------------------------------------------------------------------------------------------------

Here we can see that, for some genes, there is a high percentage of reads which are not aligned, or aligned to a gene different from the one the read was simulated from.

| index            | gene\_from |  perc\_not\_aligned|  perc\_not\_aligned\_to\_original|  perc\_aligned\_to\_original\_multimap|  perc\_aligned\_to\_original\_uniquely|
|:-----------------|:-----------|-------------------:|---------------------------------:|--------------------------------------:|--------------------------------------:|
| HLA-personalized | HLA-A      |                0.00|                              0.01|                                   1.22|                                  98.77|
| Reference        | HLA-A      |                4.69|                              4.77|                                  10.00|                                  80.55|
| HLA-personalized | HLA-B      |                0.00|                              0.00|                                   0.36|                                  99.64|
| Reference        | HLA-B      |                9.19|                             12.57|                                   4.06|                                  74.18|
| HLA-personalized | HLA-C      |                0.00|                              0.00|                                   0.01|                                  99.99|
| Reference        | HLA-C      |                4.09|                             13.82|                                   2.83|                                  79.27|
| HLA-personalized | HLA-DPB1   |                0.00|                              0.00|                                   0.00|                                 100.00|
| Reference        | HLA-DPB1   |                2.26|                              0.01|                                   0.00|                                  97.74|
| HLA-personalized | HLA-DQA1   |                0.00|                              0.00|                                   0.62|                                  99.37|
| Reference        | HLA-DQA1   |               17.96|                             12.57|                                   9.80|                                  59.67|
| HLA-personalized | HLA-DQB1   |                0.00|                              0.00|                                   0.01|                                  99.98|
| Reference        | HLA-DQB1   |               17.98|                              9.18|                                   2.37|                                  70.47|
| HLA-personalized | HLA-DRB1   |                0.00|                              0.01|                                   0.43|                                  99.56|
| Reference        | HLA-DRB1   |               12.06|                             22.01|                                   4.04|                                  61.89|

To further understand the read "loss" or "gain" we examined, for each HLA gene, the percentage of alignments stratified by mapped or source gene, respectively.

So, the size of each point in the plots below represents the percentage of alignments involving that gene pair.

### Percentage of aligments involving reads simulated from HLA in x axis but mapping to other gene

<img src="./plots/alignments_to_diff_gene.png" width="2400" />

### Percentage of alignments involving reads not simulated from HLA in x axis but still mapping to them

<img src="./plots/alignments_from_diff_gene.png" width="2400" />

Comparisons between indices and aligners
========================================

kallisto vs STAR-Salmon; HLA-diversity index
--------------------------------------------

### TPM

<img src="./plots/kallisto_vs_star_TPM.png" width="2000" />

### PCA-corrected expression

<img src="./plots/kallisto_vs_star_10pc.png" width="2000" />

kallisto vs STAR-Salmon; Reference chromosomes only
---------------------------------------------------

### TPM

<img src="./plots/kallisto_vs_star_PRI_TPM.png" width="2000" />

### PCA-corrected expression

<img src="./plots/kallisto_vs_star_PRI_10pc.png" width="2000" />

kallisto; HLA-diversity vs Reference chromosomes only
-----------------------------------------------------

### TPM

<img src="./plots/kallisto_imgt_vs_PRI_TPM.png" width="2000" />

### PCA-corrected expression

<img src="./plots/kallisto_imgt_vs_PRI_10pc.png" width="2000" />

STAR-Salmon; HLA-diversity vs Reference chromosomes only
--------------------------------------------------------

### TPM

<img src="./plots/star_imgt_vs_PRI_TPM.png" width="2000" />

### PCA-corrected expression

<img src="./plots/star_imgt_vs_PRI_10pc.png" width="2000" />
