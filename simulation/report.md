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
| HLA-A      |          2.3216|        20.3804|      0.0004|     0.7588|
| HLA-B      |          2.1864|        38.2640|      0.0002|     3.5494|
| HLA-C      |          2.2980|        28.8020|      0.0004|     0.8900|
| HLA-DPB1   |          1.2432|        10.4550|      0.0006|     0.6730|
| HLA-DQA1   |          1.1666|        20.8500|      0.0004|     4.7542|
| HLA-DQB1   |          1.3754|        22.5052|      0.0026|     5.2698|
| HLA-DRB1   |          1.9876|        66.8854|      0.0002|     7.8766|

Percentage of simulated reads from each HLA gene that aligned to a different reference:

| gene\_read |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:-----------|---------------:|--------------:|-----------:|----------:|
| HLA-A      |          1.8524|         5.1044|      1.3774|     3.6446|
| HLA-B      |          0.4026|        15.0642|      0.3470|    28.4666|
| HLA-C      |          0.0740|         9.8402|      0.0082|     9.8112|
| HLA-DPB1   |          0.0004|         0.0082|      0.0000|     0.0016|
| HLA-DQA1   |          0.8952|         5.9270|      0.5840|     6.1890|
| HLA-DQB1   |          0.1184|         5.7800|      0.0260|     9.4782|
| HLA-DRB1   |          0.9258|        24.4876|      0.4192|    20.5300|

The plot below shows where the reads simulated from a HLA gene (x-axis) are mapped when they aligned to a different reference.

The dot diameters represent the average percentage of mappings across all individuals.

![](./plots/diff_refs_alignments.png)

Percentage of simulated reads gained by each HLA gene (reads simulated from other references)

| gene\_ref |  kallisto\_imgt|  kallisto\_pri|  star\_imgt|  star\_pri|
|:----------|---------------:|--------------:|-----------:|----------:|
| HLA-A     |          0.0312|         2.5858|      0.0192|     1.6834|
| HLA-B     |          0.0368|         6.5590|      0.0120|     9.0920|
| HLA-C     |          0.1148|        16.8928|      0.0140|    24.5266|
| HLA-DPB1  |          0.0000|         0.0012|      0.0000|     0.0006|
| HLA-DQA1  |          0.1082|         0.3302|      0.0696|     0.1696|
| HLA-DQB1  |          0.0004|         0.0380|      0.0000|     0.0010|
| HLA-DRB1  |          0.1258|         9.3788|      0.0718|    10.4316|

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
