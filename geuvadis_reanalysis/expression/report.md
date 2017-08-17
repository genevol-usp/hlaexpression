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

kallisto vs GEM-based hlaTX
---------------------------

![](./plots/kallisto_vs_hlatx.png)

kallisto vs STAR-Salmon
-----------------------

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

![](./plots/kallisto_imgt_vs_chr.png)

STAR: HLA diversity vs reference chromosomes only
-------------------------------------------------

![](./plots/star_imgt_vs_chr.png)

### HLA diversity vs reference chromosomes only: R-square by quartile of difference to the reference allele

| locus    |  quartile|  rsq.kallisto.pca|  rsq.star.pca|  rsq.kallisto.tpm|  rsq.star.tpm|
|:---------|---------:|-----------------:|-------------:|-----------------:|-------------:|
| HLA-A    |         1|         0.3337326|     0.8649847|         0.7525339|     0.9426476|
| HLA-A    |         4|         0.4016592|     0.7513025|         0.8360531|     0.8901427|
| HLA-B    |         1|         0.2073240|     0.7647267|         0.6354513|     0.8529853|
| HLA-B    |         4|         0.2673982|     0.6830582|         0.3002050|     0.5370284|
| HLA-C    |         1|         0.1315891|     0.6168387|         0.5237036|     0.8817988|
| HLA-C    |         4|         0.3363081|     0.7742544|         0.5045220|     0.8496013|
| HLA-DQA1 |         1|         0.2448666|     0.5637763|         0.5845406|     0.8613144|
| HLA-DQA1 |         4|         0.0633651|     0.1525523|         0.3133389|     0.5799763|
| HLA-DQB1 |         1|         0.7892851|     0.8330204|         0.8274186|     0.7967978|
| HLA-DQB1 |         4|         0.7024806|     0.6933773|         0.7607583|     0.7217278|
| HLA-DRB1 |         1|         0.3883507|     0.6649790|         0.7357500|     0.8588556|
| HLA-DRB1 |         4|         0.4192566|     0.6262905|         0.5938372|     0.7528766|

**..................................................................**

**NOTE: All analyses below were carried out using kallisto-IMGT data**

**.................................................................**

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
