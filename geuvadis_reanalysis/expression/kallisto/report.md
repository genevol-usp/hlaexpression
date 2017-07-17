Report
================

Library Sizes and pseudoaligned reads
=====================================

![](./plots/library_sizes.png)

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

\*Allele calls were compared at the maximum resolution possible at each pair

| locus |  accuracy|
|:------|---------:|
| A     |      0.97|
| B     |      0.97|
| C     |      0.95|
| DQB1  |      0.95|
| DRB1  |      0.98|

Expression estimates
====================

GEM-based vs. kallisto
----------------------

![](./plots/kallisto_vs_hlatx.png)

Comparison to Geuvadis' published data
--------------------------------------

![](./plots/kallisto_vs_geuvadis.png)

Comparison of different indices
-------------------------------

### HLA diversity vs Reference chromosomes

![](./plots/kallisto_imgt_vs_chr.png)

### HLA diversity vs Reference chromosomes + Alternate haplotypes + Patches + Scaffolds

![](./plots/kallisto_imgt_vs_all.png)

### Distribution of TPM values

![](./plots/tpm_distributions.png)

ASE
===

ASE by number of genotype inference errors
------------------------------------------

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
