Report
================

Notes:
======

Index nomenclature:

-   pri: Primary assembly of the reference genome, containing the referece chromosomes and scaffolds (no alternate haplotypes)

-   imgt: the previous index supplemented with IMGT references

kallisto
========

Genotyping
----------

| locus    |  accuracy (%)|
|:---------|-------------:|
| HLA-A    |           100|
| HLA-B    |           100|
| HLA-C    |           100|
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
| HLA-A      |            2.30|          17.90|           0|       0.77|
| HLA-B      |            2.16|          37.31|           0|       3.53|
| HLA-C      |            2.31|          23.71|           0|       0.93|
| HLA-DQA1   |            1.13|          17.93|           0|       4.61|
| HLA-DQB1   |            1.37|          17.75|           0|       4.30|
| HLA-DRB1   |            1.98|          66.64|           0|       7.48|

Percentage of simulated reads from each HLA gene that aligned to a different reference:

| gene\_read |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:-----------|---------------:|--------------:|-----------:|----------:|
| HLA-A      |            1.98|           4.11|        1.38|       3.54|
| HLA-B      |            0.37|          13.93|        0.32|      28.51|
| HLA-C      |            0.07|           6.10|        0.01|       9.64|
| HLA-DQA1   |            0.84|           4.60|        0.56|       6.05|
| HLA-DQB1   |            0.12|           3.74|        0.02|       7.58|
| HLA-DRB1   |            0.94|          19.94|        0.45|      21.28|

The plot belows shows where the reads simulated from a HLA gene (x-axis) are mapped in alignments to a different reference (`n` in the entire dataset).

**This is just an initial exploration**

**The plot is hard to interpret, both because of the y unit and amount of colors**

![](./plots/diff_refs_alignments.png)

Percentage of simulated reads gained by each HLA gene (reads simulated from other references)

| gene\_ref |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:----------|---------------:|--------------:|-----------:|----------:|
| HLA-A     |            0.03|           1.95|        0.02|       1.56|
| HLA-B     |            0.04|           5.39|        0.01|       8.82|
| HLA-C     |            0.11|          12.78|        0.01|      25.21|
| HLA-DQA1  |            0.11|           0.00|        0.08|       0.00|
| HLA-DQB1  |            0.00|           0.00|        0.00|       0.00|
| HLA-DRB1  |            0.13|           0.00|        0.08|       0.00|

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
