Report
================

Library Sizes and pseudoaligned reads (kallisto)
================================================

![](./plots/library_sizes.png)

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

\*Allele calls were compared at the maximum resolution possible at each pair

| locus |  kallisto|  star|
|:------|---------:|-----:|
| A     |      0.97|  0.97|
| B     |      0.98|  0.98|
| C     |      0.95|  0.96|
| DQB1  |      0.95|  0.95|
| DRB1  |      0.98|  0.98|

Expression estimates
====================

kallisto vs STAR-Salmon
-----------------------

### TPM

![](./plots/kallisto_vs_star_TPM.png)

### PCA-corrected

![](./plots/kallisto_vs_star.png)

kallisto vs Geuvadis (Published)
--------------------------------

![](./plots/kallisto_vs_geuvadis.png)

STAR-Salmon vs Geuvadis (Published)
-----------------------------------

![](./plots/star_vs_geuvadis.png)

kallisto vs Geuvadis (new quantifications)
------------------------------------------

![](./plots/kallisto_vs_geuvadis_new.png)

STAR-Salmon vs Geuvadis (new quantifications)
---------------------------------------------

![](./plots/star_vs_geuvadis_new.png)

kallisto: HLA diversity vs reference chromosomes only
-----------------------------------------------------

### TPM

![](./plots/kallisto_imgt_vs_chr_TPM.png)

### PCA-corrected

![](./plots/kallisto_imgt_vs_chr.png)

STAR: HLA diversity vs reference chromosomes only
-------------------------------------------------

### TPM

![](./plots/star_imgt_vs_chr_TPM.png)

### PCA-corrected

![](./plots/star_imgt_vs_chr.png)

### HLA diversity vs reference chromosomes only: R-square by quartile of difference to the reference allele

| locus    |  quartile|  rsq.kallisto.tpm|  rsq.kallisto.pca|  rsq.star.tpm|  rsq.star.pca|
|:---------|---------:|-----------------:|-----------------:|-------------:|-------------:|
| HLA-A    |         1|              0.75|              0.33|          0.94|          0.86|
| HLA-A    |         4|              0.84|              0.40|          0.89|          0.75|
| HLA-B    |         1|              0.64|              0.21|          0.85|          0.76|
| HLA-B    |         4|              0.30|              0.27|          0.54|          0.68|
| HLA-C    |         1|              0.52|              0.13|          0.88|          0.62|
| HLA-C    |         4|              0.50|              0.34|          0.85|          0.77|
| HLA-DQA1 |         1|              0.58|              0.24|          0.86|          0.56|
| HLA-DQA1 |         4|              0.31|              0.06|          0.58|          0.15|
| HLA-DQB1 |         1|              0.83|              0.79|          0.80|          0.83|
| HLA-DQB1 |         4|              0.76|              0.70|          0.72|          0.69|
| HLA-DRB1 |         1|              0.74|              0.39|          0.86|          0.66|
| HLA-DRB1 |         4|              0.59|              0.42|          0.75|          0.63|

Distribution of TPM values
--------------------------

![](./plots/tpm_distributions.png)

ASE
===

ASE by number of genotyping errors
----------------------------------

\*Each point represents a heterozygous genotype in the intersect with Gourraud data.

\*There are more points with extreme ASE associated with genotyping errors because I'm not applying anymore a threshold of expression between 2nd allele/1st allele after the second round of the pipeline.

![](./plots/ase.png)

ASE distribution
----------------

![](./plots/ase_histogram.png)

Correlation of expression
=========================

Correlation decreases with the increase in the number of PEER factors/PCs
-------------------------------------------------------------------------

![](./plots/correlation_decrease.png)

\*Expression data in the plots below correspond to TPM values corrected by 10 PCs

Among the HLA genes
-------------------

![](./plots/hlacorrelations.png)

Between Class II genes and CIITA
--------------------------------

![](./plots/trans_activ_corrs.png)

Between pairs of HLA genes on the same vs on different haplotypes
-----------------------------------------------------------------

### HLA-A vs HLA-B

![](./plots/a_vs_b.png)

### HLA-A vs HLA-C

![](./plots/a_vs_c.png)

### HLA-B vs HLA-C

![](./plots/b_vs_c.png)

### HLA-DQA1 vs HLA-DQB1

![](./plots/dqa_vs_dqb.png)

### HLA-DQA1 vs HLA-DRB1

![](./plots/dqa_vs_drb.png)
