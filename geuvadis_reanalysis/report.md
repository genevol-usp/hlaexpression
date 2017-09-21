Report
================

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

| locus |   star|  kallisto|
|:------|------:|---------:|
| A     |  97.55|     97.06|
| B     |  98.04|     97.71|
| C     |  96.73|     96.08|
| DQB1  |  98.37|     98.53|
| DRB1  |  99.51|     99.51|

Expression estimates
====================

![](./expression/plots/expression_boxplot.png)

kallisto vs STAR-Salmon
-----------------------

### TPM

![](./expression/plots/star_vs_kallisto_TPM.png)

### PCA-corrected

![](./expression/plots/star_vs_kallisto_PCA.png)

HLA diversity vs reference chromosomes only
-------------------------------------------

### TPM

#### STAR

![](./expression/plots/star_imgt_vs_pri_TPM.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_TPM.png)

### PCA-corrected

#### STAR

![](./expression/plots/star_imgt_vs_pri_PCA.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_PCA.png)

Distribution of TPM values
--------------------------

![](./expression/plots/tpm_distributions.png)

ASE
---

### ASE by number of genotyping errors

\*Each point represents a heterozygous genotype in the intersect with Gourraud data.

![](./expression/plots/ase.png)

### ASE distribution

![](./expression/plots/ase_histogram.png)

Correlation of expression
-------------------------

![](./expression/plots/correlation_decrease.png)

### Among the HLA genes

![](./expression/plots/hlacorrelations.png)

### Between Class II genes and CIITA

![](./expression/plots/trans_activ_corrs.png)

### Between pairs of HLA genes on the same vs on different haplotypes

#### HLA-A vs HLA-B

![](./expression/plots/a_vs_b.png)

#### HLA-A vs HLA-C

![](./expression/plots/a_vs_c.png)

#### HLA-B vs HLA-C

![](./expression/plots/b_vs_c.png)

#### HLA-DQA1 vs HLA-DQB1

![](./expression/plots/dqa_vs_dqb.png)

#### HLA-DQA1 vs HLA-DRB1

![](./expression/plots/dqa_vs_drb.png)

eQTLs
=====

**All analyses were carried out using European individuals only**

PCA of genotypes
----------------

![](./qtls/plots/genotype_pca.png)

Number of eGenes according to index
-----------------------------------

![](./qtls/plots/n_of_egenes.png)

Distribution of eQTLs around the TSS
------------------------------------

### IMGT index

![](./qtls/plots/qtls_landscape_imgt.png)

### Reference transcriptome

![](./qtls/plots/qtls_landscape_pri.png)

HLA lineages
------------

![](./qtls/plots/lineage_and_effects.png)

Comparison with Geuvadis
------------------------

### RTC

Some of our eQTLs seem to mark the same biological signal as the eQTLs found by Geuvadis (RTC &gt; 0.9)

**coming soon**

### LD between our best SNP (STAR) and Geuvadis best SNP

**coming soon**

### slope and p-value

**coming soon**

Association to GWAS traits
==========================

**coming soon**

Trans-eQTLs
===========

-   Approximate pass as described on QTLtools website

**coming soon**
