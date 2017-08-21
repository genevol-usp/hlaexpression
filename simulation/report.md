Report
================

kallisto
========

Genotyping
----------

| locus |  accuracy (%)|
|:------|-------------:|
| A     |           100|
| B     |           100|
| C     |           100|
| DQA1  |           100|
| DQB1  |           100|
| DRB1  |           100|

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

| gene\_read |  kallisto|  star|
|:-----------|---------:|-----:|
| HLA-A      |      2.44|     0|
| HLA-B      |      2.35|     0|
| HLA-C      |      2.47|     0|
| HLA-DQA1   |      1.23|     0|
| HLA-DQB1   |      1.41|     0|
| HLA-DRB1   |      2.11|     0|

Percentage of simulated reads from each HLA gene that aligned to a different reference:

| gene\_read |  kallisto|  star|
|:-----------|---------:|-----:|
| HLA-A      |      2.32|  1.66|
| HLA-B      |      0.32|  0.20|
| HLA-C      |      0.08|  0.01|
| HLA-DQA1   |      0.66|  0.40|
| HLA-DQB1   |      0.05|  0.01|
| HLA-DRB1   |      0.93|  0.41|

Percentage of simulated reads gained by each HLA gene (reads simulated from other references)

| gene\_ref |  kallisto|  star|
|:----------|---------:|-----:|
| HLA-A     |      0.03|  0.02|
| HLA-B     |      0.03|  0.01|
| HLA-C     |      0.11|  0.02|
| HLA-DQA1  |      0.04|  0.03|
| HLA-DQB1  |      0.00|  0.00|
| HLA-DRB1  |      0.10|  0.05|

Comparisons between indices and aligners
========================================

kallisto vs STAR-Salmon; HLA-diversity index
--------------------------------------------

### Counts

![](./plots/kallisto_vs_star_counts.png)

### TPM

![](./plots/kallisto_vs_star_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_vs_star_10pc.png)

kallisto vs STAR-Salmon; Reference chromosomes only
---------------------------------------------------

### Counts

![](./plots/kallisto_vs_star_CHR_counts.png)

### TPM

![](./plots/kallisto_vs_star_CHR_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_vs_star_CHR_10pc.png)

kallisto; HLA-diversity vs Reference chromosomes only
-----------------------------------------------------

### Counts

![](./plots/kallisto_imgt_vs_chr_counts.png)

### TPM

![](./plots/kallisto_imgt_vs_chr_TPM.png)

### PCA-corrected expression

![](./plots/kallisto_imgt_vs_chr_10pc.png)

STAR-Salmon; HLA-diversity vs Reference chromosomes only
--------------------------------------------------------

### Counts

![](./plots/star_imgt_vs_chr_counts.png)

### TPM

![](./plots/star_imgt_vs_chr_TPM.png)

### PCA-corrected expression

![](./plots/star_imgt_vs_chr_10pc.png)

Varying the mismatch threshold in the alignment to REF chromosomes
==================================================================

TPM
---

![](./plots/star_mismatch_rates_TPM.png)

PCA-corrected expression
------------------------

![](./plots/star_mismatch_rates_PCA.png)

Mean Absolute Relative Difference
=================================

Relative difference is =

-   0 if x<sub>i</sub> = y<sub>i</sub> = 0

-   |x<sub>i</sub> - y<sub>i</sub>| / (x<sub>i</sub> + y<sub>i</sub>) otherwise

, where x<sub>i</sub> and y<sub>i</sub> are the true and estimated gene counts respectively. Then I take the mean over all genes, obtaining a value of MRD for each sample.

![](./plots/mrd.png)
