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

![](./plots/kallisto_prop_mapped.png)

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

![](./plots/star_prop_mapped.png)

Quality assessment
==================

Percentage of simulated reads not aligned:

| gene\_read |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:-----------|---------------:|--------------:|-----------:|----------:|
| HLA-A      |            2.32|          20.38|           0|       0.76|
| HLA-B      |            2.19|          38.26|           0|       3.55|
| HLA-C      |            2.30|          28.80|           0|       0.89|
| HLA-DPB1   |            1.24|          10.46|           0|       0.67|
| HLA-DQA1   |            1.17|          20.85|           0|       4.75|
| HLA-DQB1   |            1.38|          22.51|           0|       5.27|
| HLA-DRB1   |            1.99|          66.89|           0|       7.88|

Percentage of simulated reads from each HLA gene that aligned to a different reference:

| gene\_read |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:-----------|---------------:|--------------:|-----------:|----------:|
| HLA-A      |            1.85|           5.10|        1.38|       3.64|
| HLA-B      |            0.40|          15.06|        0.35|      28.47|
| HLA-C      |            0.07|           9.84|        0.01|       9.81|
| HLA-DPB1   |            0.00|           0.01|        0.00|       0.00|
| HLA-DQA1   |            0.90|           5.93|        0.58|       6.19|
| HLA-DQB1   |            0.12|           5.78|        0.03|       9.48|
| HLA-DRB1   |            0.93|          24.49|        0.42|      20.53|

The plot below shows where the reads simulated from a HLA gene (x-axis) are mapped when they aligned to a different reference.

The dot diameters represent the average percentage of mappings across all individuals.

![](./plots/diff_refs_alignments.png)

Percentage of simulated reads gained by each HLA gene (reads simulated from other references)

| gene\_ref |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:----------|---------------:|--------------:|-----------:|----------:|
| HLA-A     |            0.03|           2.59|        0.02|       1.68|
| HLA-B     |            0.04|           6.56|        0.01|       9.09|
| HLA-C     |            0.11|          16.89|        0.01|      24.53|
| HLA-DPB1  |            0.00|           0.00|        0.00|       0.00|
| HLA-DQA1  |            0.11|           0.33|        0.07|       0.17|
| HLA-DQB1  |            0.00|           0.04|        0.00|       0.00|
| HLA-DRB1  |            0.13|           9.38|        0.07|      10.43|

Comparisons between indices and aligners
========================================

kallisto vs STAR-Salmon; HLA-diversity index
--------------------------------------------

### TPM

![](./plots/kallisto_vs_star_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_vs_star_10pc.png)

kallisto vs STAR-Salmon; Reference chromosomes only
---------------------------------------------------

### TPM

![](./plots/kallisto_vs_star_PRI_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_vs_star_PRI_10pc.png)

kallisto; HLA-diversity vs Reference chromosomes only
-----------------------------------------------------

### TPM

![](./plots/kallisto_imgt_vs_PRI_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_imgt_vs_PRI_10pc.png)

STAR-Salmon; HLA-diversity vs Reference chromosomes only
--------------------------------------------------------

### TPM

![](./plots/star_imgt_vs_PRI_TPM.png)

### PCA-corrected expression

![](./plots/star_imgt_vs_PRI_10pc.png)
